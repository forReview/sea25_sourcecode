#include "external/CLI11.hpp"
#include "merge_grams.h"
#include "grammar_algorithms.h"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads=1;
    size_t n_chunks{};
    off_t chunk_size{};
    float i_frac=0.0;
    bool ver=false;
    bool rand_acc=false;
    bool skip_simp=false;
    bool skip_rl=false;
    bool part=false;
    bool check_gram=false;
    size_t seed=0;

    std::string p_file;
    std::vector<std::string> grammars_to_merge;
    std::vector<std::string> position_list;

    std::string coord_file;
    std::string version= "v1.0.1 alpha";
    //bool se_par_rounds = false;
    //size_t str{};
    //off_t start{};
    //off_t end{};
    std::vector<std::string> ra_positions;
};

struct CellWidthValidator : public CLI::Validator {
    CellWidthValidator() {
        name_ = "ValidCellWidth";
        func_ = [](const std::string &str) {
            bool valid = (str=="1") || (str=="2") || (str=="4") || (str=="8");
            if(!valid)
                return std::string(str+" is not a valid number of bytes for a native integer type");
            else
                return std::string();
        };
    }
};
const static CellWidthValidator ValidCellWidth;

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

/*
std::function<std::string(std::string)> unit_scale = [](std::string str) -> std::string {

    std::vector<std::string> pats_mb = {"mb", "MB", "mbs", "MBs", "m", "M"};
    std::vector<std::string> pats_kb = {"kb", "KB", "kbs", "KBs", "k", "K"};
    std::vector<std::string> pats_gb = {"gb", "GB", "gbs", "GBs", "g", "G"};
    std::string ext;

    size_t size, scale=0, val;

    if(ends_with(str, pats_mb, ext)){
        scale = 1048576;
    }else if(ends_with(str, pats_gb, ext)){
        scale = 1073741824;
    } else if(ends_with(str, pats_kb, ext)){
        scale = 1024;
    }

    if(!ext.empty()){
        str = str.substr(0, str.length() - ext.length());
        std::stringstream stream(str);
        stream >> val;
        size = size_t(val*scale);
        str = std::to_string(size);
    }
    return str;
};*/

static void parse_app(CLI::App& app, struct arguments& args){

    auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(27);
    app.formatter(fmt);
    app.add_flag("-v,--version", args.ver, "Print the software version and exit");

    CLI::App* comp = app.add_subcommand("cmp", "Compress text");

    //compression
    comp->add_option("TEXT", args.input_file, "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    comp->add_option("-o,--output-file", args.output_file, "Output file")->type_name("");
    comp->add_option("-t,--threads", args.n_threads, "Maximum number of parsing threads")->default_val(1);

    //comp->add_option("-s,--seed", args.seed, "Seed to generate the grammar (def. 0)");
    comp->add_flag("-q,--skip-simp", args.skip_simp, "Do not simplify the grammar");
    comp->add_flag("-e,--skip-rl", args.skip_rl, "Do not perform run-length compression");
    comp->add_flag("-r,--random-support", args.rand_acc, "Augment the grammar with random access support");
    comp->add_flag("-g,--check-gram", args.check_gram, "Check that the grammar was compressed correctly");
    //comp->add_flag("-p,--partial", args.part, "Build a partial grammar representation");

    //comp->add_option("-c,--text-chunks", args.n_chunks, "Number of text chunks in memory during the parsing (def. n_threads+1)")->default_val(0);
    comp->add_option("-f,--fraction", args.i_frac, "The parsing threads will try to use at most this input fraction");
    comp->add_option("-c,--chunk-size", args.chunk_size, "Size in bytes of each text chunk (def. min(TEXT_SIZE*0.005, 200MB))")->default_val(0);

    //metadata
    CLI::App* meta = app.add_subcommand("met", "Get the metadata of a grammar");
    meta->add_option("GRAM", args.input_file, "Input grammar in LCG format")->check(CLI::ExistingFile)->required();

    //merge
    //CLI::App* merge = app.add_subcommand("mrg", "Merge grammars");
    //merge->add_option("GRAM LIST", args.grammars_to_merge, "Grammars to be merged")->check(CLI::ExistingFile)->required();
    //merge->add_option("-o,--output-file", args.output_file, "Output grammar")->required(true);
    //merge->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")->check(CLI::ExistingDirectory)->default_val("/tmp");

    //access
    CLI::App* access = app.add_subcommand("acc", "Random access");
    access->add_option("GRAM", args.input_file, "Input grammar")->check(CLI::ExistingFile)->required();
    access->add_option("STR:START-END", args.ra_positions, "Coordinates to be accessed");

    //edit
    //CLI::App* edt = app.add_subcommand("edt", "Edit grammar-compressed text");
    //edt->add_option("GRAM", args.input_file, "Grammar to be edited")->check(CLI::ExistingFile)->required();
    //edt->add_option("STR:START-END", args.ra_positions, "Coordinates to be removed");
    //edt->add_option("-o,--output-file", args.output_file, "Output grammar");
    app.require_subcommand(1,1);
}

void comp_int(std::string& input_file, arguments& args) {
    if(args.skip_rl){
        if(args.rand_acc){
            build_gram<lc_gram_t<false, false, true>>(input_file, args.output_file, args.n_threads,
                                                       args.chunk_size, args.i_frac, args.skip_simp,
                                                       args.part, args.check_gram);
        }else{

            build_gram<lc_gram_t<false, false, false>>(input_file, args.output_file, args.n_threads,
                                                       args.chunk_size, args.i_frac, args.skip_simp,
                                                       args.part, args.check_gram);
        }
    }else{
        if(args.rand_acc){
            build_gram<lc_gram_t<false, true, true>>(input_file, args.output_file, args.n_threads,
                                                     args.chunk_size, args.i_frac, args.skip_simp,
                                                     args.part, args.check_gram);
        } else {
            build_gram<lc_gram_t<false, true, false>>(input_file, args.output_file, args.n_threads,
                                                      args.chunk_size, args.i_frac, args.skip_simp,
                                                      args.part, args.check_gram);
        }
    }
}

void access_int(std::string& input_file, std::vector<str_coord_type>& query_coords, bool has_rl_rules, bool has_cg_rules, bool has_rand_access) {
    assert(has_rand_access);
    if(has_cg_rules){
        if(has_rl_rules){
            lc_gram_t<true, true, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.start<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }else{
            lc_gram_t<true, false, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }
    }else{
        if(has_rl_rules){
            lc_gram_t<false, true, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }else{
            lc_gram_t<false, false, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }
    }
}


std::vector<str_coord_type> parse_query_coords(std::vector<std::string>& str_queries){

    size_t str;
    off_t start, end;
    std::vector<str_coord_type> query_coords;
    query_coords.reserve(str_queries.size());

    for(auto const &coord: str_queries){
        std::vector<std::string> tmp = split(coord, ':');
        if(tmp.size()!=2){
            std::cout<<"Coordinate error: query \""<<coord<<"\" is il-formed"<<std::endl;
            exit(1);
        }
        std::vector<std::string> tmp2 = split(tmp[1], '-');
        if(tmp2.size()!=2){
            std::cout<<"Coordinate error: query \""<<coord<<"\" is il-formed"<<std::endl;
            exit(1);
        }

        try{
            str = stoi(tmp[0]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp[0] <<"\" in \""<<coord<<"\" is not a valid string ID\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp[0] <<"\" in \""<<coord<<"\" is not a valid string ID\n";
            exit(1);
        }

        try{
            start = stoi(tmp2[0]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp2[0] <<"\" in \""<<coord<<"\" is not a valid start\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp2[0] <<"\" in \""<<coord<<"\" is not a valid start\n";
            exit(1);
        }

        try{
            end = stoi(tmp2[1]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp2[1] <<"\" in \""<<coord<<"\" is not a valid end\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp2[1] <<"\" in \""<<coord<<"\" is not a valid end\n";
            exit(1);
        }
        query_coords.emplace_back(str, start, end);
    }
    return query_coords;
}

int main(int argc, char** argv) {

    arguments args;
    CLI::App app("Locally-consistent grammar compression");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(args.ver){
        std::cout<<args.version<<std::endl;
        exit(0);
    }

    if(app.got_subcommand("cmp")) {
        std::cout << "\nInput file: " << args.input_file << " ("<<report_space(file_size(args.input_file))<<")"<<std::endl;
        if (args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
        args.output_file = std::filesystem::path(args.output_file).replace_extension(".lcg");
        std::string input_collection = args.input_file;
        comp_int(input_collection, args);
    } else if(app.got_subcommand("met")){
        print_metadata(args.input_file);
    } else if (app.got_subcommand("acc")){
        std::vector<str_coord_type> query_coords = parse_query_coords(args.ra_positions);
        bool has_rl_rules, has_cg_rules, has_rand_access;
        std::tie(has_rl_rules, has_cg_rules, has_rand_access) = read_grammar_flags(args.input_file);
        access_int(args.input_file, query_coords, has_rl_rules, has_cg_rules, has_rand_access);
    } else {
        std::cout<<" Unknown command "<<std::endl;
        return 1;
    }
    return 0;
}
