#define RANGELESS_FN_ENABLE_PARALLEL 1
#include "fn.hpp"



/////////////////////////////////////////////////////////////////////////////
namespace gx
{

class thread_pool
{
public:
    thread_pool(uint32_t max_threads_ = std::thread::hardware_concurrency(),
                uint32_t capacity_    = std::thread::hardware_concurrency() * 16)

        : m_max_threads{ max_threads_ }
        , m_queue{ capacity_ }
    {
        // if either is 0, will execute in this thread
    }

    /////////////////////////////////////////////////////////////////////////

    uint32_t max_threads() const noexcept
    {
        return m_max_threads;
    }

    uint32_t capacity() const noexcept
    {
        return (uint32_t)m_queue.capacity();
    }

    /// Number of pending + executing tasks.
    size_t size() const noexcept
    {
        // note: this is not just the size of the pending queue:
        // even if that is 0, there may be some tasks still
        // executing, so we keep track of size ourselves.
        return m_tasks_count;
    }

    /// If true and no other thread is enqueueing, indicates that the
    /// pending queue is empty and all tasks have been completed.
    bool empty() const noexcept
    {
        return m_tasks_count == 0;
    }

    /// Obtains the lock; waits until the jobs complete, and resets the pool.
    void flush()
    {
        std::lock_guard<std::mutex> g{ m_mutex };

        // Note that submitting a single, or one-per-worker empty
        // is not sufficient, because there may be a long-running
        // job(s), and the barrier task(s) are quickly processed by
        // a single worker and their futures will become ready
        // while the bona-fide jobs are still running.

        // Submit empty stopper tasks.
        for(size_t i = 0; i < m_workers.size(); ++i) {
            m_queue.push( task_t{} );
        }

        // Wait for each worker to finish.
        for(std::future<void>& fut : m_workers) {
            fut.wait();
        }

        // Rethrow first exception, if any.
        // (if that happens, other tasks are not running because
        // we waited on them in the above loop)
        for(std::future<void>& fut : m_workers) {
            fut.get();
        }

        m_workers.clear();

        //assert(m_queue.empty());
        assert(m_tasks_count == 0);
    }


    /////////////////////////////////////////////////////////////////////////
    /*! \brief Enqueue a job (nullary invokable)
     *
     * NB: this method can be called from multiple threads.
     * Jobs will be executed in order of submission.
     * Blocks if the jobs-queue is at capacity.
     * Returned future is non-blocking in destructor.
    */
    template<typename F>
    auto operator()(F fn) -> std::future< decltype( fn() )>
    {
        using ret_t = decltype( fn() );
        auto task = s_package_task<ret_t>(std::move(fn));
        auto ret = task.get_future();

        if(m_max_threads == 0 || m_queue.capacity() == 0) {
            // complete the task in this thread.
            assert(m_tasks_count == 0);
            task();
            return ret;
        }

        std::lock_guard<std::mutex> g{ m_mutex };

        assert(m_max_threads > 0);
        assert(m_workers.size() <= m_max_threads);

        const bool spawn_new_worker =
            m_workers.empty()                   ? true  // need at least one thread
          : m_workers.size() >= m_max_threads   ? false // have enough threads
          : m_workers.size() > m_tasks_count    ? false // have a waiting thread to do it
          :                                       true; // need another thread

        ++m_tasks_count; //corresponding decrement is in the main_loop

        if(spawn_new_worker) {
            m_workers.emplace_back(
                std::async(
                    std::launch::async,
                    &thread_pool::x_main_loop,
                    this));
        }

        m_queue.push( s_package_task<void>( std::move(task) ) );

        return ret;
    }


    /////////////////////////////////////////////////////////////////////////
    // \brief Calls flush().
    ~thread_pool()
    {
        flush();
        assert(m_workers.empty());
        //assert(m_queue.empty());
        assert(m_tasks_count == 0);
    }

    /////////////////////////////////////////////////////////////////////////
private:
    using  task_t = std::packaged_task<void()>;
    using queue_t = rangeless::mt::synchronized_queue<task_t>;

                   const uint32_t m_max_threads;
    std::deque<std::future<void>> m_workers;
                          queue_t m_queue;
                       std::mutex m_mutex; // protects m_workers
               std::atomic_size_t m_tasks_count = { 0 };

    /////////////////////////////////////////////////////////////////////////
private:

    template<typename ret_t, typename F>
    static auto s_package_task(F fn) -> std::packaged_task<ret_t()>
    {
#ifdef _MSC_VER
        // MSVC bug: can't create packaged-task if F is move-only;
        // as a work-around, wrap fn in a shared-pointer and capture
        // it in a closure.
        auto fn_ptr = std::make_shared<F>( std::move(fn) );
        return std::packaged_task<ret_t()>{ [fn_ptr] { return (*fn_ptr)(); } };
#else
        return std::packaged_task<ret_t()>{ std::move(fn) };
#endif
    }

    /////////////////////////////////////////////////////////////////////////
    // worker-thread loop: process tasks until a stopper task
    // (pushed by destructor) is encountered.
    void x_main_loop() noexcept
    {
        for(auto task = m_queue.pop();
                 task.valid();
                 task = m_queue.pop())
        {
            task();
            --m_tasks_count;
        }
    }

}; // pool

}
