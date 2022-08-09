#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>

#include "strips.h"
#include "utils.h"

#include <dfa.h>

using namespace std;

vector<string> g_default_regularizers{ };
vector<string> g_default_symmetries{ "ordering-arities-for-atoms", "ordering-action-arguments", "ordering-objects-in-layers" };
vector<string> g_default_encoding_options{ "applicable-actions-tuples", "features-to-actions-map-strict-lex-ordering" };

inline string color(const string &control_sequence, const string &str, bool use_colors = true) {
    return use_colors ? control_sequence + str + Utils::normal() : str;
}
inline string cyan(const string &str, bool use_colors = true) {
    return color(Utils::cyan(), str, use_colors);
}
inline string magenta(const string &str, bool use_colors = true) {
    return color(Utils::magenta(), str, use_colors);
}

void read_model(ostream &os, const string &filename, map<string, vector<vector<string> > > &model, bool disable_colors, bool debug) {
    os << cyan(string("reading file '") + filename + "' ...", !disable_colors) << flush;
    ifstream ifs(filename.c_str());
    if( !ifs.fail() ) {
        for( string line; getline(ifs, line); ) {
            size_t index = line.find_first_of(')');
            assert(index != string::npos);
            string key = line.substr(0, 1 + index);

            // extract tuple
            string rest = line.substr(1 + index);
            if( rest[0] == '(' ) {
                if( debug ) os << "record: key=" << key << ", tuple=" << rest << ", values={";
                vector<vector<string> > &records = model[key];
                vector<string> tuple;
                char *cp = strdup(&rest[1]);
                for( char *p = strtok(cp, ",)"); p != nullptr; p = strtok(nullptr, ",)") ) {
                    tuple.push_back(p);
                    if( debug ) os << tuple.back() << ",";
                }
                if( debug ) os << "}" << endl;
                records.emplace_back(move(tuple));
            }
        }
        os << cyan(" done!", !disable_colors) << endl;
        ifs.close();
    } else {
        os << Utils::error() << "opening file '" << filename << "'" << endl;
        exit(-1);
    }
}

template<typename T>
const DFA::DFA<T>* read_dfa(const string &dfa_filename, vector<string> &comments, bool verbose, bool use_colors = true) {
    if( verbose ) cout << cyan(string("reading '") + dfa_filename + "' ... ", use_colors) << flush;
    const DFA::DFA<T>* dfa = nullptr;
    ifstream ifs(dfa_filename.c_str());
    if( !ifs.fail() ) {
        dfa = DFA::DFA<T>::read_dump(ifs, comments);
        ifs.close();
        if( verbose ) cout << cyan(" done!", use_colors) << endl;
    } else {
        throw runtime_error(Utils::error() + "opening DFA file '" + dfa_filename + "'");
    }
    return dfa;
}

template<typename T>
const DFA::DFA<T>* read_dfa(const string &dfa_filename, bool verbose, bool use_colors = true) {
    vector<string> comments;
    return read_dfa<T>(dfa_filename, comments, verbose, use_colors);
}

string call(int argc, const char **argv) {
    string call_str;
    for( int i = 0; i < argc; ++i ) {
        call_str += argv[i];
        if( 1 + i < argc ) call_str += " ";
    }
    return call_str;
}

void usage(ostream &os, const string &name) {
    string pad = string(string("usage: ").length() + name.length(), ' ');
    os << endl
       << "usage: " << name << " --dump-ts-dot [--disable-colors] [--help] [--output <filename>] <dfa>" << endl
       << endl
       << "       " << name << " [--amo-encoding { Quad | Log | Heule }]" << endl
       << pad + " [--debug]" << endl
       << pad + " [--decode-full]" << endl
       << pad + " [--decode-meta-layer]" << endl
       << pad + " [--decode-model]" << endl
       << pad + " [--disable-colors]" << endl
       << pad + " [--disable-symmetry-handling]" << endl
       << pad + " [--encoding { <option> | -<option> | clear }]" << endl
       << pad + " [--help]" << endl
       << pad + " [--input <filename>]" << endl
       << pad + " [--max-number-ground-actions <layer> <n>]" << endl
       << pad + " [--output <filename>]" << endl
       << pad + " [--partial-assignment <assignment>]" << endl
       << pad + " [--regularizer { <option> | -<option> | clear }]" << endl
       << pad + " [--symmetries { <option> | -<option> | clear }]" << endl
       << pad + " [--verify-meta-layer <meta-layer>]" << endl
       << pad + " [--weak-lex-orderings]" << endl
       << pad + " [--weighted-sat]" << endl
       << pad + " <num-action-schemas> <max-arity-first-action> ... <max-arity-last-action>" << endl
       << pad + " <num-atom-schemas> <max-arity-first-atom> ... <max-arity-last-atom>" << endl
       << pad + " <num-meta-features> <num-static-unary-predicates> <num-static-binary-predicates>" << endl
       << pad + " [<num-features> <num-objects> <LB-on-sum-features-per-state> <UB-on-sum-features-per-state> <prefix>]*" << endl
       << endl
       << "where (all required!)" << endl
       << "    <num-action-schemas> is number of action schemas in template" << endl
       << "    <num-meta-features> is number of meta-features for actions (partitioned between actions)" << endl
       << "    <num-meta-objects> is number of meta-objects for actions (define arguments and are shared between actions)" << endl
       << "    <num-atom-schemas> is number of available atom schemas (meta-features are mapped into atoms)" << endl
       << "    <max-arity> is bound for maximum arity of atoms" << endl
       << "    <num-static-unary-predicates> is bound for number of static unary predicates" << endl
       << "    <num-static-binary-predicates> is bound for number of static binary predicates" << endl
       << endl
       << "where (required for each layer!)" << endl
       << "    <num-features> is number of boolean features in layer" << endl
       << "    <num-objects> is number of objects in layer" << endl
       << "    <LB-on-sum-features-per-state> is lower bound on sum of true features across all states (0 means disabled)" << endl
       << "    <UB-on-sum-features-per-state> is upper bound on sum of true features across all states (0 means disabled)" << endl
       << "    <prefix> is prefix for layer (used to get .dfa and other info for layer)" << endl
       << endl
       << "For the options," << endl
       << "    --amo-encoding to set the type of encoding for AMO boolean constraints" << endl
       << "    --debug to turn on debugging" << endl
       << "    --decode-full to *fully* decode model" << endl
       << "    --decode-meta-layer to extract meta-layer in model" << endl
       << "    --decode-model to decode model" << endl
       << "    --disable-colors to disable colors in output" << endl
       << "    --disable-symmetry-handling do not generate formulas for reducing symmetries" << endl
       << "    --dump-ts-dot to dump transition system DFA in dot format" << endl;

    using Theory = StructuralLearning::Theory;

    // encoding options
    os << "    --encoding <option> to insert into current set, -<option> to remove from current set, 'clear' to erase current set (it may be used more than once)" << endl;
    os << "        " << Utils::bold() << "available:" << Utils::normal();
    vector<string> available_encoding_options = Theory::get_available_encoding_options();
    if( available_encoding_options.empty() ) os << " <none>";
    for( int i = 0; i < int(available_encoding_options.size()); ++i )
        os << " " << available_encoding_options.at(i);
    os << endl;
    os << "        " << Utils::bold() << "default:" << Utils::normal();
    if( g_default_encoding_options.empty() ) os << " <none>";
    for( int i = 0; i < int(g_default_encoding_options.size()); ++i )
        os << " " << g_default_encoding_options.at(i);
    os << endl;

    // continue with other options
    os << "    --help to print this help and exit" << endl
       << "    --input to set filename of input file to <filename>" << endl
       << "    --max-number-ground-actions to bound number of ground actions in layer <layer> to <n>" << endl
       << "    --output to set filename of output file to <filename>" << endl
       << "    --partial-assignment to extend theory with literals in <assignment>" << endl;

    // options for regularization
    os << "    --regularizer <option> to insert into current set, -<option> to remove from current set, 'clear' to erase current set (it may be used more than once)" << endl;
    os << "        " << Utils::bold() << "available:" << Utils::normal();
    vector<string> available_regularizers = Theory::get_available_regularizers();
    if( available_regularizers.empty() ) os << " <none>";
    for( int i = 0; i < int(available_regularizers.size()); ++i )
        os << " " << available_regularizers.at(i);
    os << endl;
    os << "        " << Utils::bold() << "default:" << Utils::normal();
    if( g_default_regularizers.empty() ) os << " <none>";
    for( int i = 0; i < int(g_default_regularizers.size()); ++i )
        os << " " << g_default_regularizers.at(i);
    os << endl;

    // options for handling symmetries
    os << "    --symmetries <option> to insert into current set, -<option> to remove from current set, 'clear' to erase current set (it may be used more than once)" << endl;
    os << "        " << Utils::bold() << "available:" << Utils::normal();
    vector<string> available_symmetries = Theory::get_available_symmetries();
    if( available_symmetries.empty() ) os << " <none>";
    for( int i = 0; i < int(available_symmetries.size()); ++i )
        os << " " << available_symmetries.at(i);
    os << endl;
    os << "        " << Utils::bold() << "default:" << Utils::normal();
    if( g_default_symmetries.empty() ) os << " <none>";
    for( int i = 0; i < int(g_default_symmetries.size()); ++i )
        os << " " << g_default_symmetries.at(i);
    os << endl;

    // continue with other options
    os << "    --verify-meta-layer to verify meta-layer in <meta-layer>" << endl
       << "    --weak-lex-orderigs to enforce weak (non-strict) lexicographic orderings (rather than strong or strict)" << endl
       << "    --weighted-sat to generate partial weighted Max-SAT theory" << endl
       << endl
       ;
}

void insufficient_arguments(ostream &os, const string &name) {
    os << Utils::error() << "insufficient arguments" << endl;
    usage(os, name);
    exit(0);
}

void decode(const string &request,
            const StructuralLearning::Theory &theory,
            const map<string, string> &output_filename_map,
            const vector<string> &output_filenames,
            bool disable_colors) {
    // open output file (if any)
    ostream *os = &cout;
    if( (output_filename_map.count(request) > 0) || !output_filenames.empty() ) {
        string output_filename;
        if( output_filename_map.count(request) > 0 )
            output_filename = output_filename_map.at(request);
        else
            output_filename = output_filenames.back();
        cout << cyan(string("writing file '") + output_filename + "' ...", !disable_colors) << flush;
        os = new ofstream(output_filename.c_str());
    }

    // do decoding
    if( request == "decode-full" )
        theory.decode_model_full(*os);
    else if( request == "decode-meta-layer" )
        theory.decode_meta_layer(*os);
    else if( request == "decode-model" )
        theory.decode_model(*os);

    // close output file (if any)
    if( os != &cout ) {
        static_cast<ofstream*>(os)->close();
        delete os;
        cout << cyan(" done!", !disable_colors) << endl;
    }
}

int main(int argc, const char **argv) {
    using Theory = StructuralLearning::Theory;
    using MetaLayerParameters = StructuralLearning::MetaLayerParameters;
    using LayerParameters = StructuralLearning::LayerParameters;
    using Options = StructuralLearning::Options;

    // print call
    const string call_str(call(argc, argv));
    cout << "call: " << call_str << endl;

    // read executable name
    string executable_name(*argv);

    // options
    SAT::Theory::amo_encoding_t opt_amo_encoding = SAT::Theory::amo_encoding_t::Heule;
    bool opt_debug = false;
    bool opt_decode_full = false;
    bool opt_decode_meta_layer = false;
    bool opt_decode_model = false;
    bool opt_disable_colors = false;
    bool opt_disable_symmetry_handling = false;
    bool opt_dump_ts_dot = false;
    set<string> opt_encoding;
    string opt_input_filename;
    map<int, int> opt_max_number_ground_actions;
    vector<string> opt_output_filenames;
    map<string, string> opt_output_filename_map;
    vector<string> opt_partial_assignment_filenames;
    vector<pair<int, int> > opt_random_transitions;
    set<string> opt_regularizers;
    set<string> opt_symmetries;
    bool opt_verify_meta_layer = false;
    bool opt_weak_lex_orderings = false;
    bool opt_weighted_sat = false;
    string opt_input_meta_layer;

    // add default regularizers, symmetries and encoding options
    Theory::add_regularizers(opt_regularizers, g_default_regularizers);
    Theory::add_symmetries(opt_symmetries, g_default_symmetries);
    Theory::add_encoding_options(opt_encoding, g_default_encoding_options);

    // read options
    for( ++argv, --argc; (argc > 0) && (**argv == '-'); ++argv, --argc ) {
        if( string(*argv) == "--amo-encoding" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            string amo_encoding = *argv;
            if( amo_encoding == "Quad" ) {
                opt_amo_encoding = SAT::Theory::amo_encoding_t::Quad;
            } else if( amo_encoding == "Log" ) {
                opt_amo_encoding = SAT::Theory::amo_encoding_t::Log;
            } else if( amo_encoding == "Heule" ) {
                opt_amo_encoding = SAT::Theory::amo_encoding_t::Heule;
            } else {
                cout << Utils::error() << "unrecognized AMO encoding '" << *argv << "'" << endl;
                usage(cout, executable_name);
                exit(0);
            }
        } else if( string(*argv) == "--debug" ) {
            opt_debug = true;
        } else if( string(*argv) == "--decode-full" ) {
            opt_decode_full = true;
            if( !opt_output_filenames.empty() ) {
                opt_output_filename_map.emplace("decode-full", opt_output_filenames.back());
                opt_output_filenames.pop_back();
            }
        } else if( string(*argv) == "--decode-meta-layer" ) {
            opt_decode_meta_layer = true;
            if( !opt_output_filenames.empty() ) {
                opt_output_filename_map.emplace("decode-meta-layer", opt_output_filenames.back());
                opt_output_filenames.pop_back();
            }
        } else if( string(*argv) == "--decode-model" ) {
            opt_decode_model = true;
            if( !opt_output_filenames.empty() ) {
                opt_output_filename_map.emplace("decode-model", opt_output_filenames.back());
                opt_output_filenames.pop_back();
            }
        } else if( string(*argv) == "--disable-colors" ) {
            opt_disable_colors = true;
        } else if( string(*argv) == "--disable-symmetry-handling" ) {
            opt_disable_symmetry_handling = true;
        } else if( string(*argv) == "--dump-ts-dot" ) {
            opt_dump_ts_dot = true;
            if( !opt_output_filenames.empty() ) {
                opt_output_filename_map.emplace("dump-ts-dot", opt_output_filenames.back());
                opt_output_filenames.pop_back();
            }
        } else if( string(*argv) == "--encoding" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            string option = *argv;
            if( option == "clear" )
                Theory::clear_all_encoding_options(opt_encoding);
            else if( option.at(0) == '-' )
                Theory::clear_encoding_option(opt_encoding, option.substr(1));
            else
                Theory::add_encoding_option(opt_encoding, option);
        } else if( string(*argv) == "--help" ) {
            usage(cout, executable_name);
            exit(0);
        } else if( string(*argv) == "--input" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            opt_input_filename = *argv;
        } else if( string(*argv) == "--max-number-ground-actions" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            int layer = atoi(*argv);
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            int n = atoi(*argv);
            if( !opt_max_number_ground_actions.emplace(layer, n).second ) {
                cout << Utils::warning() << "bound for layer " << layer
                     << " already specified; ignoring new value " << n
                     << "..."
                     << endl;
            }
        } else if( string(*argv) == "--output" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            opt_output_filenames.emplace_back(*argv);
        } else if( string(*argv) == "--partial-assignment" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            opt_partial_assignment_filenames.emplace_back(*argv);
        } else if( string(*argv) == "--regularizer" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            string option = *argv;
            if( option == "clear" )
                Theory::clear_all_regularizers(opt_regularizers);
            else if( option.at(0) == '-' )
                Theory::clear_regularizer(opt_regularizers, option.substr(1));
            else
                Theory::add_regularizer(opt_regularizers, option);
        } else if( string(*argv) == "--symmetries" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            string option = *argv;
            if( option == "clear" )
                Theory::clear_all_symmetries(opt_symmetries);
            else if( option.at(0) == '-' )
                Theory::clear_symmetry(opt_symmetries, option.substr(1));
            else
                Theory::add_symmetry(opt_symmetries, option);
        } else if( string(*argv) == "--random-transitions" ) {
            if( argc < 3 ) insufficient_arguments(cout, executable_name);
            int layer = atoi(argv[1]);
            int n = atoi(argv[2]);
            opt_random_transitions.emplace_back(layer, n);
            argv += 2;
            argc -= 2;
        } else if( string(*argv) == "--verify-meta-layer" ) {
            ++argv;
            --argc;
            if( argc == 0 ) insufficient_arguments(cout, executable_name);
            opt_verify_meta_layer = true;
            opt_input_meta_layer = *argv;
            Theory::add_encoding_option(opt_encoding, "features-to-actions-map-none");
        } else if( string(*argv) == "--weak-lex-orderings" ) {
            opt_weak_lex_orderings = true;
        } else if( string(*argv) == "--weighted-sat" ) {
            opt_weighted_sat = true;
        } else {
            cout << Utils::error() << "unrecognized option '" << *argv << "'" << endl;
            usage(cout, executable_name);
            exit(0);
        }
    }

    // either make state-space dfa or SAT-related stuff
    if( opt_dump_ts_dot ) {
        // read arguments
        if( argc == 0 ) insufficient_arguments(cout, executable_name);
        string opt_dfa_filename = *argv++;
        --argc;

        // read dfa
        const DFA::DFA<string> *dfa = read_dfa<string>(opt_dfa_filename, true, !opt_disable_colors);

        // open output file (if any)
        ostream *os = &cout;
        if( (opt_output_filename_map.count("dump-ts-dot") > 0) || !opt_output_filenames.empty() ) {
            string output_filename;
            if( opt_output_filename_map.count("dump-ts-dot") > 0 )
                output_filename = opt_output_filename_map.at("dump-ts-dot");
            else
                output_filename = opt_output_filenames.back();
            cout << cyan(string("writing file '") + output_filename + "' ...", !opt_disable_colors) << flush;
            os = new ofstream(output_filename.c_str());
        }

        // dump dfa
        dfa->dump_dot(*os, !opt_disable_colors);

        // close output file (if any)
        if( os != &cout ) {
            static_cast<ofstream*>(os)->close();
            delete os;
            cout << cyan(" done!", !opt_disable_colors) << endl;
        }

        // remove dfa
        delete dfa;
    } else {
        // read arguments
        if( argc == 0 ) insufficient_arguments(cout, executable_name);
        int opt_number_actions = atoi(*argv++);
        vector<int> opt_arity_for_actions(opt_number_actions);
        if( --argc < opt_number_actions ) insufficient_arguments(cout, executable_name);
        for( int i = 0; i < opt_number_actions; ++i, --argc )
            opt_arity_for_actions[i] = atoi(*argv++);

        if( argc == 0 ) insufficient_arguments(cout, executable_name);
        int opt_number_atoms = atoi(*argv++);
        vector<int> opt_arity_for_atoms(opt_number_atoms);
        if( --argc < opt_number_atoms ) insufficient_arguments(cout, executable_name);
        for( int i = 0; i < opt_number_atoms; ++i, --argc )
            opt_arity_for_atoms[i] = atoi(*argv++);

        if( argc < 3 ) insufficient_arguments(cout, executable_name);
        int opt_number_meta_features = atoi(argv[0]);
        int opt_number_static_unary_predicates = atoi(argv[1]);
        int opt_number_static_binary_predicates = atoi(argv[2]);

        // parameters for meta layer
        MetaLayerParameters parameters_for_meta_layer(opt_arity_for_actions,
                                                      opt_arity_for_atoms,
                                                      opt_number_meta_features,
                                                      opt_number_static_unary_predicates,
                                                      opt_number_static_binary_predicates);

        // read arguments for meta layers
        vector<pair<vector<int>, string> > opt_layers;
        for( argc -= 3, argv += 3; argc != 0; argc -= 5, argv += 5 )
            opt_layers.emplace_back(vector<int>{ atoi(argv[0]), atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) }, argv[4]);

        // input filenames
        vector<LayerParameters> parameters_for_layers;
        for( int i = 0; i < int(opt_layers.size()); ++i ) {
            int num_features = opt_layers[i].first.at(0);
            int num_objects = opt_layers[i].first.at(1);
            int feature_sum_per_state_lower = opt_layers[i].first.at(2);
            int feature_sum_per_state_upper = opt_layers[i].first.at(3);
            string dfa_filename = opt_layers[i].second;
            if( (dfa_filename.size() < 4) || (dfa_filename.substr(dfa_filename.size() - 4) != ".dfa") )
                dfa_filename += ".dfa";

            // read data and create space
            const DFA::DFA<string> *dfa = read_dfa<string>(dfa_filename, true, !opt_disable_colors);
            parameters_for_layers.emplace_back(num_features, num_objects, feature_sum_per_state_lower, feature_sum_per_state_upper, dfa);
        }

        // verification of meta-layer
        if( opt_verify_meta_layer ) {
            if( !opt_partial_assignment_filenames.empty() ) {
                cout << Utils::error() << "cannot use --verify-meta-layer and --partial-assignment simultaneously; stop" << endl;
                exit(0);
            }
            opt_partial_assignment_filenames.push_back(opt_input_meta_layer);
        }

        // if verifying meta-layer, read arities from given partial assignment
        if( opt_verify_meta_layer ) {
            string partial_assignment_filename = opt_partial_assignment_filenames.front();
            cout << cyan(string("reading file '") + partial_assignment_filename + "' ...", !opt_disable_colors) << flush;
            ifstream ifs(partial_assignment_filename.c_str());
            if( !ifs.fail() ) {
                for( string line; getline(ifs, line); ) {
                    if( line.substr(0, 10) == "arity(p,i)" ) {
                        istringstream iss(line.substr(10));
                        char left, comma, right;
                        int p, i;
                        iss >> left >> p >> comma >> i >> right;
                        assert((left == '(') && (comma == ',') && (right == ')'));
                        assert((0 <= p) && (p < int(parameters_for_meta_layer.atoms_.size())));
                        assert((0 <= i) && (i < 1 + parameters_for_meta_layer.atoms_.at(p)));
                        parameters_for_meta_layer.atoms_[p] = i;
                    }
                }
                cout << cyan(" done!", !opt_disable_colors) << endl;
            } else {
                cout << Utils::error() << "opening file '" << partial_assignment_filename << "'" << endl;
                exit(-1);
            }
        }

        // set options for SAT theory
        Options options;
        options.amo_encoding_ = opt_amo_encoding;
        options.debug_ = opt_debug;
        options.decode_ = opt_decode_full || opt_decode_meta_layer || opt_decode_model;
        options.disable_colors_ = opt_disable_colors;
        options.disable_symmetry_handling_ = opt_disable_symmetry_handling;
        options.encoding_ = opt_encoding;
        options.partial_assignment_ = !opt_partial_assignment_filenames.empty();
        options.regularizers_ = opt_regularizers;
        options.symmetries_ = opt_symmetries;
        options.verify_meta_layer_ = opt_verify_meta_layer;
        options.weak_lex_orderings_ = opt_weak_lex_orderings;
        options.weighted_sat_ = opt_weighted_sat;
        options.max_number_ground_actions_ = opt_max_number_ground_actions;
        options.random_transitions_ = opt_random_transitions;
        options.normalize();

        // create SAT theory
        float start_time_build_variables = Utils::read_time_in_seconds();
        Theory *theory = new Theory(parameters_for_meta_layer, parameters_for_layers, options);
        float elapsed_time_build_variables = Utils::read_time_in_seconds() - start_time_build_variables;

        // decode
        if( opt_decode_full || opt_decode_meta_layer || opt_decode_model ) {
            // perform void construction to terminate generation of (auxiliary) variables
            theory->build_theory();

            // report theory size
            cout << "#variables=" << theory->num_variables() << endl;
            cout << "#implications=" << theory->num_implications() << endl;
            cout << "#soft-implications=" << theory->num_soft_implications() << endl;

            // read model
            string model_filename = opt_input_filename.empty() ? "model.cnf" : opt_input_filename;
            cout << cyan(string("reading file '") + model_filename + "' ...", !opt_disable_colors) << flush;
            ifstream ifs(model_filename.c_str());
            if( !ifs.fail() ) {
                // read model w/ appropriate reader; to decide, extract header and re-open file
                string header;
                ifs >> header;
                ifs.close();
                ifs.open(model_filename.c_str());
                assert(!ifs.fail());

                if( (header == "SAT") || (header == "UNSAT") )
                    theory->read_minisat_output(ifs);
                else if( isdigit(*header.c_str()) || (*header.c_str() == '-') )
                    theory->read_glucose_output(ifs);
                else
                    theory->read_other_output(ifs);

                ifs.close();
                cout << cyan(" done!", !opt_disable_colors) << endl;

                // do decoding
                if( opt_decode_full )
                    decode("decode-full", *theory, opt_output_filename_map, opt_output_filenames, opt_disable_colors);
                if( opt_decode_meta_layer )
                    decode("decode-meta-layer", *theory, opt_output_filename_map, opt_output_filenames, opt_disable_colors);
                if( opt_decode_model )
                    decode("decode-model", *theory, opt_output_filename_map, opt_output_filenames, opt_disable_colors);
            } else {
                cout << Utils::error() << "opening file '" << model_filename << "'" << endl;
            }
        } else {
            // open output stream (tunnel)
            string theory_filename = string("strips") + (opt_weighted_sat ? ".wcnf" : ".cnf");
            if( !opt_output_filenames.empty() ) theory_filename = opt_output_filenames.back();
            cout << cyan(string("set tunnel to '") + theory_filename + "'", !opt_disable_colors) << endl;
            ofstream ofs(theory_filename.c_str());
            theory->set_tunnel(&ofs, opt_weighted_sat);

            // output header stub if weighted SAT
            if( opt_weighted_sat ) {
                ofs << "c this header requires fix!" << endl;
                ofs << theory->header(opt_weighted_sat) << endl;
            }

            // build theory (generate output)
            float start_time_build_theory = Utils::read_time_in_seconds();
            theory->add_comment(string("call: ") + call_str);
            theory->build_theory();
            float elapsed_time_build_theory = Utils::read_time_in_seconds() - start_time_build_theory;

            // read partial assignments (if any)
            for( int i = 0; i < int(opt_partial_assignment_filenames.size()); ++i ) {
                cout << cyan(string("reading file '") + opt_partial_assignment_filenames[i] + "' ...", !opt_disable_colors) << flush;
                ifstream ifs(opt_partial_assignment_filenames[i].c_str());
                if( !ifs.fail() ) {
                    theory->add_comment("partial assignment");
                    pair<int, int> p = theory->read_assignment(ifs, true);
                    ifs.close();
                    cout << cyan(string(" done! (#lines=") + to_string(p.first) + ", #units=" + to_string(p.second) + ")", !opt_disable_colors) << endl;
                } else {
                    cout << Utils::warning() << "can't open file '" << opt_partial_assignment_filenames[i] << "'" << endl;
                }
            }

            // report theory size and timing statistics
            cout << "#variables=" << theory->num_variables() << endl;
            cout << "#implications=" << theory->num_implications() << endl;
            cout << "#soft-implications=" << theory->num_soft_implications() << endl;
            cout << "time-for-building-variables=" << elapsed_time_build_variables << endl;
            cout << "time-for-building-theory=" << elapsed_time_build_theory << endl;

            // close output stream (tunnel)
            if( !opt_weighted_sat ) ofs << theory->header(opt_weighted_sat) << endl;
            ofs.close();
            theory->set_tunnel(nullptr);
        }

        // free resources
        delete theory;
        for( int i = 0; i < int(parameters_for_layers.size()); ++i )
            delete parameters_for_layers[i].dfa_;
    }
    return 0;
}

