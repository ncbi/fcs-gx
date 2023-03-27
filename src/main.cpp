// See PUBLIC DOMAIN NOTICE at the bottom.

#include "util.hpp"
#include "serial_util.hpp"
#include "tmasker.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ext/CLI11.hpp"
#pragma GCC diagnostic pop

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>

#if not defined(GIT_REVISION)
#    define GIT_REVISION "unknown"
#endif
extern const char* const g_git_revision;
const char* const g_git_revision = GIT_REVISION;


// 'path' is a path to .gxi or a directory containing a unique .gxi file.
// postcondition: path is a path to .gxi; else error-out; no-op if path is empty.
static std::string resolve_path_to_gxi(const std::string& path)
{
    using namespace gx;

    if (path.empty() || str::endswith(path, ".gxi")) {
        return path;
    }

    auto dirp = opendir(path.c_str());

    if (!dirp) {
        std::cerr << "Error: Expected a path to a direcotry containing a single gx-database, or a path to the *.gxi file:  " << path << "\n";
        std::exit(1);
    }

    std::string result_path = "";
    size_t gxi_count = 0;
    while (auto dp = readdir(dirp))
        if (str::endswith(dp->d_name, ".gxi"))
    {
        ++gxi_count;
        result_path = path + "/" + dp->d_name;
    }

    if (gxi_count == 0) {
        std::cerr << "Error: Did not find a unique *.gxi file in " << path 
                  << " - it does not appear to be a gx-database dir.\n";
        std::exit(1);
    } else if (gxi_count > 1) {
        std::cerr << "Error: Found multile *.gxi files in " << path
                  << "\nPlease put databases in separate dirs, or specify a full path to a .gxi file.\n";
        std::exit(1);
    }

    (void)closedir(dirp);
    errno = 0;
    return result_path;
}


static int run(int argc, char* argv[])
{
    using namespace gx;

    ser::prepare_standard_streams();

    CLI::App app{"GX aligner and taxonomy classifier tool.\n\nhttps://github.com/ncbi/fcs/wiki/FCS-GX\n"};

    std::string seq_id2tax_id_path;
    std::string hardmask_locs_path;
    std::string softmask_locs_path;
    std::string action_report_path;
    std::string contam_fasta_out_path;
    std::string fasta_for_repeatdb_path;
    std::string taxa_path;
    std::string prot_minhash_id;
    std::string db_path   = "db/all.gxi";
    std::string inp_path  = "stdin";
    std::string out_path  = "stdout";
    std::string asserted_div = "unknown";
    size_t action_report_min_seq_len = 200;

    auto add_inp_arg = [&inp_path](CLI::App& cmd_, const std::string& descr)
    {
        // cmd_ - trailing underscore business because false-positive shadow warning.
        cmd_.add_option("-i,--input", inp_path)
            ->type_name("file")
            ->capture_default_str()
            ->description(descr + "\n");
    };

    auto add_out_arg = [&out_path](CLI::App& cmd_, const std::string& descr = "")
    {
        cmd_.add_option("-o,--output", out_path)
            ->type_name("file")
            ->capture_default_str()
            ->description(descr + "\n");

        // --help first by default; move it at the bottom.
        cmd_.set_help_flag("-h,--help")
            ->description("Print this help message and exit.\n");
    };

    auto add_db_opt = [&db_path](CLI::App& cmd_)
    {
        cmd_.add_option("--gx-db", db_path)
            ->type_name("file")
            ->capture_default_str()
            ->description("/dev/shm/path/to/gxdb/ produced by gx make-db.\nShould be placed in RAM-disk.\n");
    };

    auto add_repeats_opt = [&fasta_for_repeatdb_path](auto& cmd_)
    {
        cmd_.add_option("--repeats-basis-fa", fasta_for_repeatdb_path)
            ->type_name("file")
            ->description("Fasta to gather transposon-repeats statistics (normally same input as for --input).\n");
    };

    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("make-db", "Create GX-database from a set of subject sequences.");

        add_inp_arg(cmd, "fasta");

        cmd.add_option("--seq_id-tax_id", seq_id2tax_id_path)
            ->required(true)
            ->type_name("file")
            ->description("2-column TSV with the following header:\n"
                          GX_TSV_HEADER__SEQ_ID_MAPPING "\n"
                          "#seq-id\\ttax-id\n");

        cmd.add_option("--taxa", taxa_path)
            ->required(true)
            ->type_name("file")
            ->description("5-column TSV with the following header:\n"
                          GX_TSV_HEADER__TAXA "\n"
                          "#tax-id\\tspecies\\tcommon-name\\tBLAST-div\\tdiv\n");

        cmd.add_option("--hardmask", hardmask_locs_path)
            ->type_name("file")
            ->description("3-column TSV with the following header:\n"
                          GX_TSV_HEADER__LOCS "\n"
                          "#seq-id\\tfrom1\\tto1\n");

        cmd.add_option("--softmask", softmask_locs_path)
            ->type_name("file")
            ->description("Same format as --hardmask.\n"
                          "Exclude locations from indexing, but do not hardmask the sequence.\n");

        add_out_arg(cmd, "/path/to/db.gxi          - the db-index file (must have .gxi suffix).\n"
                         "/path/to/db.gxs          - additional output: sequences file.\n"
                         "/path/to/db.seq_info.tsv - additional output: sequences metadata.\n"
                         "/path/to/db.meta.jsonl   - additional output: db-metadata.\n");
    }
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("get-fasta", "Fetch fasta from db.");

        add_inp_arg(cmd, "3-column TSV of locs, same format as --hardmask in make-db mode.\n");
        add_db_opt(cmd);
        add_out_arg(cmd, "Fasta.");
    }
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("split-fasta", "Split input fasta on N-runs of length at least 10.");

        add_inp_arg(cmd, "fasta");
        add_out_arg(cmd, "Fasta with seq-ids transformed as $seq_id~$start..$stop (1-based, inclusive).");
    }
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("align", "Align DNA-queries against a GX-database.");

        add_inp_arg(cmd, "fasta");
        add_db_opt(cmd);
        add_repeats_opt(cmd);
        add_out_arg(cmd, "Tabular report of ungapped alignments.");
    }
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("taxify", "Generate taxonomy-report for query sequences.");
        add_inp_arg(cmd, "hits-TSV emitted by `gx align`");
        add_db_opt(cmd);

        cmd.add_option("--asserted-div", asserted_div)
            ->required(false)
            ->type_name("gx-tax-div")
            ->capture_default_str()
            ->description("E.g. 'prok:firmicutes'");

        // Could call this option --hardmask to be consistent with make-db mode,
        // but a user may erroneously interpret it as pertaining to query.
        cmd.add_option("--db-exclude-locs", hardmask_locs_path)
            ->required(false)
            ->type_name("file")
            ->capture_default_str()
            ->description("Ignore db-hits at these locations. Same format as --hardmask arg in make-db mode.");

        add_out_arg(cmd, "Tabular taxonomy report.");
    }
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("clean-genome", "Apply edit-actions to the original input fasta.");
        add_inp_arg(cmd, "Original fasta.");

        cmd.add_option("--action-report", action_report_path)
            ->required(true)
            ->type_name("file")
            ->description("Action-report (*.fcs_gx_report.txt produced by `run_gx`).");

        cmd.add_option("--contam-fasta-out", contam_fasta_out_path)
            ->type_name("file")
            ->required(false)
            ->description("Fasta for EXCLUDE entries in the action-report.");

        cmd.add_option("--min-seq-len", action_report_min_seq_len)
            ->required(false)
            ->type_name("int")
            ->capture_default_str()
            ->description("Minimumm sequence length to keep.");

        add_out_arg(cmd, "Modified (cleaned) fasta.");
    } 
    // ----------------------------------------------------------------------
    {
        auto& cmd = *app.add_subcommand("find-repeats", "Find transposon locations.");
        add_inp_arg(cmd, "fasta");
        add_repeats_opt(cmd);
        add_out_arg(cmd, "Tabular report of locations.");
    }

    // "Hidden" commands below - only exposed in NCBI environment.
    // ----------------------------------------------------------------------
    if (std::getenv("NCBI")) {
        auto& cmd = *app.add_subcommand("prot-minhash-create", "Generate proteome-minhash.");
        add_inp_arg(cmd, "Proteome fasta");
        cmd.add_option("--prot-minhash-id", prot_minhash_id)->required(true);
        add_out_arg(cmd, "Minhash signature.");
    }
    // ----------------------------------------------------------------------
    if (std::getenv("NCBI")) {
        auto& cmd = *app.add_subcommand("prot-minhash-compare", "Pairwise-compare minhashes.");
        add_inp_arg(cmd, "Minhashes produced by gx prot-minhash-create command.");
        add_out_arg(cmd, "Pairwise distance report.");
    }
    // ----------------------------------------------------------------------
    {
        app.add_subcommand("show-license", "Show licensing information.");
    }

    app.require_subcommand(1);
    app.footer("\nbuild:" __DATE__ " " __TIME__ "; git:" + std::string{g_git_revision} + "\n");
    CLI11_PARSE(app, argc, argv);
    errno = 0; // clean-up after parse, just in case

    const auto& command = app.get_subcommands().front()->get_name();

    VERIFY(!inp_path.empty());
    VERIFY(!out_path.empty());

    std::istream& istr = inp_path == "stdin" ? std::cin :
        [&]() -> std::ifstream&
        {
            static auto ifstr = open_ifstream(inp_path);
            return ifstr;
        }();

    if (command == "make-db") {
        auto seq_ids_ifstr      = open_ifstream(seq_id2tax_id_path);
        auto taxa_ifstr_ptr     = open_ifstream_opt(taxa_path);
        auto hardmask_ifstr_ptr = open_ifstream_opt(hardmask_locs_path);
        auto softmask_ifstr_ptr = open_ifstream_opt(softmask_locs_path);

        VERIFY(out_path != "stdout");
        MakeDb(istr,
               seq_ids_ifstr,
               taxa_ifstr_ptr.get(),
               hardmask_ifstr_ptr.get(),
               softmask_ifstr_ptr.get(),
               out_path);
        return 0;
    }

    std::ostream& ostr = out_path == "stdout" ? std::cout :
        [&]() -> std::ofstream&
        {
            static std::ofstream ofstr{ out_path };
            return ofstr;
        }();

    db_path = resolve_path_to_gxi(db_path);

    if (command == "align") {
        open_ifstream(db_path); // just to throw with good error message if not accessible
        ProcessQueries(db_path, taxa_path, fasta_for_repeatdb_path, istr, ostr);

    } else if (command == "get-fasta") {
        open_ifstream(db_path); // just to throw with good error message if not accessible
        GetFasta(db_path, istr, ostr);

    } else if (command == "find-repeats") {
        auto fasta_istr1 = open_ifstream(fasta_for_repeatdb_path);
        CTmasker(fasta_istr1).process_fasta(istr, ostr);

    } else if (command == "prot-minhash-create") {
        ostr << prot_minhash_id << "\t" << MakeProtsetMinhash(istr) << "\n";

    } else if (command == "prot-minhash-compare") {
        PairwiseCompareMinHashes(istr, ostr);

    } else if (command == "taxify") {
        auto taxa_ifstr = open_ifstream(str::replace(db_path, ".gxi", ".taxa.tsv"));
        Taxify(istr, taxa_ifstr, hardmask_locs_path, asserted_div, db_path, out_path);

    } else if (command == "clean-genome") {
        auto action_report_ifstr = open_ifstream(action_report_path);

        auto contam_fasta_out_ofstr = 
            contam_fasta_out_path.empty() ? std::unique_ptr<std::ofstream>{}
                                          : std::make_unique<std::ofstream>(contam_fasta_out_path);

        ApplyActionReport(istr, action_report_ifstr, ostr, contam_fasta_out_ofstr.get(), action_report_min_seq_len);

    } else if (command == "split-fasta") {
        SplitFasta(istr, ostr);

    } else if (command == "show-license") {
        std::cout << GET_EMBEDDED_BLOB(______LICENSE);
        //   underscores correspond to ../../LICENSE in linker command in CMakeLists.txt
    } else  {
        VERIFY(false);
        return -1;
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    static const bool enable_main_try_catch = gx::get_env("GX_ENABLE_MAIN_TRY_CATCH", 1);

    if (enable_main_try_catch) try {
        // Do catch exceptions when we want the stack unwinding to happen,
        // which may add additional relevant info during the process.
        return run(argc, argv);

    } catch(const std::exception& e) {
        std::cerr << "Fatal error (" << typeid(e).name() << "): " << e.what() << "\n";
        return EXIT_FAILURE;
    } // catch(...) - not doing that because have nothing to report
      //            - let std::terminate to happen naturally

    // Don't catch exceptions when we want libdw + libbackward to get the full stack-trace,
    // but the stack-unwinding will not happen after std::terminate is called by c++ runtime.
    return run(argc, argv);
}

/////////////////////////////////////////////////////////////////////////////


/*  $Id$
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*/
