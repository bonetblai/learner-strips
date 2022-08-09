#ifndef STRIPS_H
#define STRIPS_H

#include <cassert>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <sat/encoder/theory.h>
#include <dfa.h>
#include "utils.h"

// strong typedefs for safer implementation
#include <boost/serialization/strong_typedef.hpp>

namespace StructuralLearning {

inline std::string str_color(const std::string &control_sequence, const std::string &str, bool use_colors = true) {
    return use_colors ? control_sequence + str + Utils::normal() : str;
}
inline std::string cyan(const std::string &str, bool use_colors = true) {
    return str_color(Utils::cyan(), str, use_colors);
}
inline std::string green(const std::string &str, bool use_colors = true) {
    return str_color(Utils::green(), str, use_colors);
}
inline std::string magenta(const std::string &str, bool use_colors = true) {
    return str_color(Utils::magenta(), str, use_colors);
}
inline std::string red(const std::string &str, bool use_colors = true) {
    return str_color(Utils::red(), str, use_colors);
}
inline std::string yellow(const std::string &str, bool use_colors = true) {
    return str_color(Utils::yellow(), str, use_colors);
}

inline std::string alt(bool use_colors = true) {
    return cyan("(alt)", use_colors);
}
inline std::string alt_default(bool use_colors = true) {
    return cyan("(alt-default)", use_colors);
}

inline std::string new_formula(bool use_colors = true) {
    return magenta("(new)", use_colors);
}
inline std::string updated(bool use_colors = true) {
    return magenta("(updated)", use_colors);
}
inline std::string subsumed(bool use_colors = true) {
    return magenta("(subsumed)", use_colors);
}
inline std::string disabled(bool use_colors = true) {
    return magenta("(disabled)", use_colors);
}

struct Options {
    SAT::Theory::amo_encoding_t amo_encoding_;
    bool debug_;
    bool decode_;
    bool disable_colors_;
    bool disable_symmetry_handling_;
    std::set<std::string> encoding_;           // options for encoding
    bool partial_assignment_;
    std::set<std::string> regularizers_;       // set of active regularizers
    std::set<std::string> symmetries_;         // options for handling symmetries
    bool verify_meta_layer_;
    bool weak_lex_orderings_;
    bool weighted_sat_;

    std::map<int, int> max_number_ground_actions_;
    std::vector<std::pair<int, int> > random_transitions_;

    Options()
      : amo_encoding_(SAT::Theory::amo_encoding_t::Heule),
        debug_(false),
        decode_(false),
        disable_colors_(false),
        disable_symmetry_handling_(false),
        partial_assignment_(false),
        verify_meta_layer_(false),
        weak_lex_orderings_(false),
        weighted_sat_(false) {
    }
    Options(const Options &opt)
      : amo_encoding_(opt.amo_encoding_),
        debug_(opt.debug_),
        decode_(opt.decode_),
        disable_colors_(opt.disable_colors_),
        disable_symmetry_handling_(opt.disable_symmetry_handling_),
        encoding_(opt.encoding_),
        partial_assignment_(opt.partial_assignment_),
        regularizers_(opt.regularizers_),
        symmetries_(opt.symmetries_),
        verify_meta_layer_(opt.verify_meta_layer_),
        weak_lex_orderings_(opt.weak_lex_orderings_),
        weighted_sat_(opt.weighted_sat_),
        max_number_ground_actions_(opt.max_number_ground_actions_),
        random_transitions_(opt.random_transitions_) {
    }
    void normalize() { }
};

// create strong typedefs (type aliases)
BOOST_STRONG_TYPEDEF(int, Layer)
BOOST_STRONG_TYPEDEF(int, Action)
BOOST_STRONG_TYPEDEF(int, MetaFeature)
BOOST_STRONG_TYPEDEF(int, Atom)
BOOST_STRONG_TYPEDEF(int, MetaObject)
BOOST_STRONG_TYPEDEF(int, Arity)
BOOST_STRONG_TYPEDEF(int, Feature)
BOOST_STRONG_TYPEDEF(int, Object)
BOOST_STRONG_TYPEDEF(int, State)
BOOST_STRONG_TYPEDEF(int, Transition)
BOOST_STRONG_TYPEDEF(int, Label)
BOOST_STRONG_TYPEDEF(int, Unary)
BOOST_STRONG_TYPEDEF(int, Binary)

struct MetaLayerParameters {
    const std::vector<int> actions_;
    const int num_actions_;
    int num_meta_objects_;
    std::vector<int> atoms_; // it may be modified during enc0
    const int num_atoms_;
    int max_arity_;
    const int num_meta_features_;
    const int num_static_unary_predicates_;
    const int num_static_binary_predicates_;
    MetaLayerParameters(const std::vector<int> &actions,
                        const std::vector<int> &atoms,
                        int num_meta_features,
                        int num_static_unary_predicates,
                        int num_static_binary_predicates)
      : actions_(actions),
        num_actions_(actions.size()),
        atoms_(atoms),
        num_atoms_(atoms.size()),
        num_meta_features_(num_meta_features),
        num_static_unary_predicates_(num_static_unary_predicates),
        num_static_binary_predicates_(num_static_binary_predicates) {
        num_meta_objects_ = 0;
        for( int i = 0; i < int(actions_.size()); ++i )
            num_meta_objects_ = std::max(num_meta_objects_, actions_[i]);
        max_arity_ = 0;
        for( int i = 0; i < int(atoms_.size()); ++i )
            max_arity_ = std::max(max_arity_, atoms_[i]);
    }
};

struct LayerParameters {
    const int num_features_;
    const int num_objects_;
    const int feature_sum_per_state_lower_;
    const int feature_sum_per_state_upper_;
    const DFA::DFA<std::string> *dfa_;
    LayerParameters(int num_features,
                    int num_objects,
                    int feature_sum_per_state_lower,
                    int feature_sum_per_state_upper,
                    const DFA::DFA<std::string> *dfa)
      : num_features_(num_features),
        num_objects_(num_objects),
        feature_sum_per_state_lower_(feature_sum_per_state_lower),
        feature_sum_per_state_upper_(feature_sum_per_state_upper),
        dfa_(dfa) {
    }
};

class Theory : public SAT::Theory {
  protected:
    const Action num_actions_;                                    // number of actions
    const std::vector<int> arity_for_actions_;                    // arity for each action CHECK: exact or upper bound?
    const MetaObject num_meta_objects_;                           // number of meta-objects
    const Atom num_atoms_;                                        // number of atom schemas
    const std::vector<int> arity_for_atoms_;                      // arity for each atom CHECK: exact or upper bound?
    const Arity max_arity_;                                       // maximum arity for atoms
    const MetaFeature num_meta_features_;                         // number of meta-features
    const Unary num_static_unary_predicates_;                     // number of static unary predicates
    const Binary num_static_binary_predicates_;                   // number of static binary predicates
    const Layer num_layers_;                                      // number of layers
    const Options options_;                                       // options

    // parameters for variables
    std::vector<Layer> p_layers_;
    std::vector<Action> p_actions_;
    std::vector<MetaFeature> p_meta_features_;
    std::vector<Atom> p_atoms_;
    std::vector<MetaObject> p_meta_objects_;
    std::vector<Arity> p_arities_;
    std::vector<Feature> p_features_;
    std::vector<Object> p_objects_;
    std::vector<State> p_states_;
    std::vector<Transition> p_transitions_;
    std::vector<Label> p_labels_;
    std::vector<Unary> p_unary_;
    std::vector<Binary> p_binary_;

    Feature num_features_;                                        // number of features across layers
    State num_states_;                                            // number of states across layers
    Transition num_transitions_;                                  // number of transitions across layers
    Label num_labels_;                                            // number of transition labels across layers

    Object max_num_objects_;                                      // max number of objects in layer
    std::vector<Object> objects_per_layer_;                       // number of objects in each layer
    std::vector<int> feature_sum_per_layer_lower_;                // (constant) sum of true features in each state per layer LB
    std::vector<int> feature_sum_per_layer_upper_;                // (constant) sum of true features in each state per layer UB

    std::vector<std::string> label_as_string_;                    // labels
    std::map<std::string, Label> label_map_;                      // labels in all ts

    std::vector<std::vector<Feature> > features_per_layer_;       // features at each layer
    std::vector<std::vector<State> > states_per_layer_;           // states at each layer
    std::vector<std::vector<Transition> > transitions_per_layer_; // transitions at each layer

    // features (each feature is unique)
    std::vector<Layer> f_layer_;

    // states (each state is unique)
    std::vector<Layer> s_layer_;

    // transitions (each transition is unique)
    std::vector<Layer> tr_layer_;
    std::vector<State> tr_src_;
    std::vector<State> tr_dst_;
    std::vector<Label> tr_label_;

    // classes (sets) of variables
    SAT::VarSet pre0, pre1;                  // [DECISION] preconditions for all fluents
    SAT::VarSet eff0, eff1;                  // [DECISION] effects for boolean fluents
    SAT::VarSet label;                       // [DECISION] action labels
    SAT::VarSet using1, using2;              // [IMPLIED] used meta-features
    SAT::VarSet arity;                       // [DECISION] arity of atoms
    SAT::VarSet atom1, atom2;                // [DECISION] mapping of meta-features to atoms, and meta-objetcs to args
    SAT::VarSet non_static0, non_static1;    // [IMPLIED] atoms must be non-static
    SAT::VarSet unary, binary;               // [DECISION] use of unary and binary static relations by actions
    SAT::VarSet args, relevant;              // [IMPLIED]

    SAT::VarSet map, mapf, phi;
    SAT::VarSet mapt, g, free;
    SAT::VarSet appl, mapeq, eq;
    SAT::VarSet Z0, Z1, X0, X1;

    SAT::VarSet ground1, ground2;
    SAT::VarSet r, s;

    SAT::VarSet W;                 // for more compact [ mapf(t,k,mk) & atom(mk,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ] ]
    SAT::VarSet gdiff;             // for encodings 1 and 2
    SAT::VarSet gdiff1, gdiff2;    // for encoding 3

    SAT::VarSet U, B;
    SAT::VarSet ord;

    // FULL (default) ENCODING OF GROUND ACTIONS:
    SAT::VarSet gtuple, G;
    SAT::VarSet appl2, violated0, violated1, pre0eq, pre1eq, eq2;
    //SAT::VarSet pos, neg;

    // size of random subset of transitions at each layer (-1 = all)
    std::vector<int> random_transitions_per_layer_;

    // formula groups
    std::vector<std::string> group_descriptions_;

  public:
    Theory(const MetaLayerParameters &parameters_for_meta_layer,
           const std::vector<LayerParameters> &parameters_for_layers,
           const Options &options)
      : SAT::Theory(options.decode_, options.amo_encoding_),
        num_actions_(parameters_for_meta_layer.actions_.size()),
        arity_for_actions_(parameters_for_meta_layer.actions_),
        num_meta_objects_(parameters_for_meta_layer.num_meta_objects_),
        num_atoms_(parameters_for_meta_layer.num_atoms_),
        arity_for_atoms_(parameters_for_meta_layer.atoms_),
        max_arity_(parameters_for_meta_layer.max_arity_),
        num_meta_features_(parameters_for_meta_layer.num_meta_features_),
        num_static_unary_predicates_(parameters_for_meta_layer.num_static_unary_predicates_),
        num_static_binary_predicates_(parameters_for_meta_layer.num_static_binary_predicates_),
        num_layers_(parameters_for_layers.size()),
        options_(options) {
        initialize(std::cout, parameters_for_layers);
        build_variables();
    }
    virtual ~Theory() { }

    static const std::vector<std::string> get_available_regularizers() {
        return std::vector<std::string>{ "exact-arities", "disjoint-meta-features" };
    }
    static void add_regularizer(std::set<std::string> &regularizers, const std::string &reg) {
        regularizers.insert(reg);
    }
    static void add_regularizers(std::set<std::string> &regularizers, const std::vector<std::string> &regs) {
        for( int i = 0; i < int(regs.size()); ++i )
            add_regularizer(regularizers, regs[i]);
    }
    static void clear_regularizer(std::set<std::string> &regularizers, const std::string &reg) {
        regularizers.erase(reg);
    }
    static void clear_all_regularizers(std::set<std::string> &regularizers) {
        regularizers.clear();
    }
    bool regularizer(const std::string &reg) const {
        return options_.regularizers_.find(reg) != options_.regularizers_.end();
    }

    static const std::vector<std::string> get_available_symmetries() {
        return std::vector<std::string>{
            "ordering-arities-for-atoms",
            "ordering-action-arguments",
            "ordering-objects-in-layers",
            "strict-ordering-objects-in-layers"
        };
    }
    static void add_symmetry(std::set<std::string> &symmetries, const std::string &symm) {
        symmetries.insert(symm);
    }
    static void add_symmetries(std::set<std::string> &symmetries, const std::vector<std::string> &symms) {
        for( int i = 0; i < int(symms.size()); ++i )
            add_symmetry(symmetries, symms[i]);
    }
    static void clear_symmetry(std::set<std::string> &symmetries, const std::string &symm) {
        symmetries.erase(symm);
    }
    static void clear_all_symmetries(std::set<std::string> &symmetries) {
        symmetries.clear();
    }
    bool symmetries(const std::string &option) const {
        return !options_.disable_symmetry_handling_ && (options_.symmetries_.find(option) != options_.symmetries_.end());
    }
    bool ordering_objects_in_layers() const {
        return symmetries("ordering-objects-in-layers") || symmetries("strict-ordering-objects-in-layers");
    }

    static const std::vector<std::string> get_available_encoding_options() {
        return std::vector<std::string>{
            "applicable-actions-tuples",
            "applicable-actions-alt",
            "features-to-actions-map-strict-lex-ordering",
            "features-to-actions-map-alt1",
            "features-to-actions-map-alt2",
            "features-to-actions-map-alt3",
            "features-to-actions-map-none"
        };
    }
    static void add_encoding_option(std::set<std::string> &options, const std::string &opt) {
        if( opt == "applicable-actions-tuples" ) {
            clear_encoding_option(options, "applicable-actions-alt");
        } else if( opt == "applicable-actions-alt" ) {
            clear_encoding_option(options, "applicable-actions-tuples");
        } else if( opt == "features-to-actions-map-strict-lex-ordering" ) {
            add_encoding_option(options, "features-to_actions-map-none");
        } else if( opt == "features-to-actions-map-alt1" ) {
            add_encoding_option(options, "features-to_actions-map-none");
        } else if( opt == "features-to-actions-map-alt2" ) {
            add_encoding_option(options, "features-to_actions-map-none");
        } else if( opt == "features-to-actions-map-alt3" ) {
            add_encoding_option(options, "features-to_actions-map-none");
        } else if( opt == "features-to-actions-map-none" ) {
            clear_encoding_option(options, "features-to_actions-map-strict-lex-ordering");
            clear_encoding_option(options, "features-to_actions-map-alt1");
            clear_encoding_option(options, "features-to_actions-map-alt2");
            clear_encoding_option(options, "features-to_actions-map-alt3");
        }
        options.insert(opt);
    }
    static void add_encoding_options(std::set<std::string> &options, const std::vector<std::string> &opts) {
        for( int i = 0; i < int(opts.size()); ++i )
            add_encoding_option(options, opts[i]);
    }
    static void clear_encoding_option(std::set<std::string> &options, const std::string &opt) {
        options.erase(opt);
    }
    static void clear_all_encoding_options(std::set<std::string> &options) {
        options.clear();
    }
    bool encoding(const std::string &option) const {
        return options_.encoding_.find(option) != options_.encoding_.end();
    }
    bool some_features_to_actions_map_encoding() const {
        return encoding("features-to-actions-map-strict-lex-ordering") ||
          encoding("features-to-actions-map-alt1") ||
          encoding("features-to-actions-map-alt2") ||
          encoding("features-to-actions-map-alt3");
    }

  protected:
    void initialize(std::ostream &os, const std::vector<LayerParameters> &parameters_for_layers) {
        os << "Theory: parameters:"
           << " #layers=" << num_layers_
           << ", #actions=" << num_actions_
           << ", #meta-objects=" << num_meta_objects_
           << ", #meta-features=" << num_meta_features_
           << ", #atoms=" << num_atoms_
           << ", max-arity=" << max_arity_
           << ", #static-unary-predicates=" << num_static_unary_predicates_
           << ", #static-binary-predicates=" << num_static_binary_predicates_
           << ", TS=[(#states,#edges,#labels)*]=["
           << std::flush;

        for( int i = 0; i < int(parameters_for_layers.size()); ++i ) {
            const DFA::DFA<std::string> &dfa = *parameters_for_layers[i].dfa_;
            os << "(" << dfa.num_states()
               << "," << dfa.num_edges()
               << "," << dfa.num_labels()
               << ")," << std::flush;
        }
        os << "]" << std::flush;

        // extract number random transitions per layer from options
        // value -1 means 'all transitions'
        random_transitions_per_layer_ = std::vector<int>(num_layers_, -1);
        for( int i = 0; i < int(options_.random_transitions_.size()); ++i ) {
            int layer = options_.random_transitions_[i].first;
            int n = options_.random_transitions_[i].second;
            random_transitions_per_layer_.at(layer) = n;
        }

        // extract/consolidate basic DFA info
        num_features_ = 0;
        num_states_ = 0;
        num_transitions_ = 0;

        // CHECK: all TS must have identical set of labels!
        // CHECK: watchout case of disjoint-meta-features: arity for actions/atoms
        // CHECK: output info during initialization
        // CHECK: default options: enc4, non-disjoint-meta-features, read-dfa
        // CHECK: encoding 0, verify meta-layer

        assert(num_layers_ == int(parameters_for_layers.size()));
        bool using_enc0 = false;
        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            const LayerParameters &parameters = parameters_for_layers[layer];
            const DFA::DFA<std::string> &dfa = *parameters.dfa_;

            // features
            int num_features_in_layer = 0;
            if( !some_features_to_actions_map_encoding() || (parameters.num_features_ == 0) ) {
                if( !using_enc0 && (layer > 0 ) )
                    throw std::runtime_error(Utils::error() + "use of enc0 must be consistent across all layers");
                using_enc0 = true;

                os << yellow("Encoding 0:", !options_.disable_colors_) << " "
                   << "bypassing given num-features in layer " << layer << " (value=" << parameters.num_features_ << ")"
                   << std::endl;

                int num_objects = parameters.num_objects_;
                for( Atom p = Atom(0); p < num_atoms_; ++p ) {
                    int arity = arity_for_atoms_[p];
                    int num_ground_atoms = 1;
                    for( int i = 0; i < arity; ++i )
                        num_ground_atoms *= num_objects;
                    num_features_in_layer += num_ground_atoms;
                }

                os << yellow("Encoding 0:", !options_.disable_colors_) << " "
                   << "using " << num_features_in_layer << " features in layer " << layer
                   << " (complete grounding of atoms)"
                   << std::endl;
            } else {
                num_features_in_layer = parameters.num_features_;
            }

            features_per_layer_.emplace_back(num_features_in_layer);
            std::iota(features_per_layer_.back().begin(), features_per_layer_.back().end(), num_features_);
            num_features_ += num_features_in_layer;
            for( int i = 0; i < num_features_in_layer; ++i )
                f_layer_.push_back(layer);

            if( options_.debug_ ) {
                os << "Features: ";
                std::copy(features_per_layer_.back().begin(),
                          features_per_layer_.back().end(),
                          std::ostream_iterator<int>(os, " "));
                os << std::endl;
            }

            // states
            int states_base = num_states_;
            states_per_layer_.emplace_back(dfa.num_states());
            std::iota(states_per_layer_.back().begin(), states_per_layer_.back().end(), num_states_);
            num_states_ += dfa.num_states();
            for( int i = 0; i < dfa.num_states(); ++i )
                s_layer_.push_back(layer);

            if( options_.debug_ ) {
                os << "States: ";
                std::copy(states_per_layer_.back().begin(),
                          states_per_layer_.back().end(),
                          std::ostream_iterator<int>(os, " "));
                os << std::endl;
            }

            // labels
            for( int i = 0; i < int(dfa.labels().size()); ++i ) {
                const std::string &label = dfa.labels()[i];
                if( label_map_.find(label) == label_map_.end() ) {
                    label_map_.emplace(label, Label(label_map_.size()));
                    label_as_string_.push_back(label);
                }
            }

            // read edges (transitions)
            std::vector<int> edge_pool(dfa.num_edges(), 0);
            std::iota(edge_pool.begin(), edge_pool.end(), 0);
            std::vector<bool> chosen_edges(dfa.num_edges(), false);

            int num_selected_edges = -1;
            int n = random_transitions_per_layer_.at(layer);
            if( (n == -1) || (dfa.num_edges() <= n) ) {
                num_selected_edges = dfa.num_edges();
                random_transitions_per_layer_.at(layer) = -1;
                chosen_edges = std::vector<bool>(dfa.num_edges(), true);
            } else {
                std::default_random_engine generator;
                num_selected_edges = n;
                while( n > 0 ) {
                    assert(!edge_pool.empty());
                    std::uniform_int_distribution<int> dist(0, edge_pool.size() - 1);
                    int index = dist(generator);
                    chosen_edges[edge_pool[index]] = true;
                    edge_pool[index] = edge_pool.back();
                    edge_pool.pop_back();
                    --n;
                }
            }

            // transitions
            transitions_per_layer_.emplace_back(num_selected_edges);
            std::iota(transitions_per_layer_.back().begin(), transitions_per_layer_.back().end(), num_transitions_);
            num_transitions_ += num_selected_edges;

            if( options_.debug_ ) {
                os << "Transitions: ";
                std::copy(transitions_per_layer_.back().begin(),
                          transitions_per_layer_.back().end(),
                          std::ostream_iterator<int>(os, " "));
                os << std::endl;
            }

            // extract transitions
            int edge_index = 0;
            int num_processed = 0;
            for( int src = 0; src < dfa.num_states(); ++src ) {
                const std::vector<std::pair<int, int> > &edges = dfa.edges(src);
                for( int i = 0; i < int(edges.size()); ++i ) {
                    if( chosen_edges.at(edge_index) ) {
                        int label_index = edges[i].first;
                        int dst = edges[i].second;
                        Label label = label_map_[dfa.get_label(label_index)];
                        tr_layer_.push_back(layer);
                        tr_src_.push_back(State(states_base + src));
                        tr_dst_.push_back(State(states_base + dst));
                        tr_label_.push_back(label);
                        ++num_processed;
                    }
                    ++edge_index;
                }
            }
            assert(num_processed == num_selected_edges);

            if( options_.debug_ ) {
                os << "dfa[" << layer << "]="
                   << "(" << dfa.num_states()
                   << "," << dfa.num_edges()
                   << ",labels={";
                std::copy(dfa.labels().begin(),
                          dfa.labels().end(),
                          std::ostream_iterator<std::string>(os, ","));
                os << "})" << std::endl;
            }
        }
        num_labels_ = label_map_.size();

        os << ", labels={";
        std::for_each(label_map_.begin(),
                      label_map_.end(),
                      [&os](const std::pair<std::string, int> &p) {
                          os << p.second << "." << p.first << "," << std::flush;
                      });
        os << "}" << std::flush;
        os << std::endl;

        // objects in layer
        max_num_objects_ = 0;
        for( Layer l = Layer(0); l < num_layers_; ++l ) {
            Object num_objects(parameters_for_layers[l].num_objects_);
            objects_per_layer_.push_back(num_objects);
            max_num_objects_ = std::max(max_num_objects_, num_objects);
        }

        // feature sums per layer
        feature_sum_per_layer_lower_ = std::vector<int>(num_layers_);
        feature_sum_per_layer_upper_ = std::vector<int>(num_layers_);
        for( Layer l = Layer(0); l < num_layers_; ++l ) {
            feature_sum_per_layer_lower_[l] = parameters_for_layers[l].feature_sum_per_state_lower_;
            feature_sum_per_layer_upper_[l] = parameters_for_layers[l].feature_sum_per_state_upper_;
        }

        // construct parameter lists
        p_layers_ = std::vector<Layer>(num_layers_);
        std::iota(p_layers_.begin(), p_layers_.end(), 0);
        if( options_.debug_ ) std::cout << "layers: " << p_layers_.size() << std::endl;

        p_atoms_ = std::vector<Atom>(num_atoms_);
        std::iota(p_atoms_.begin(), p_atoms_.end(), 0);
        p_arities_ = std::vector<Arity>(1 + max_arity_);
        std::iota(p_arities_.begin(), p_arities_.end(), 0);
        p_meta_features_ = std::vector<MetaFeature>(num_meta_features_);
        std::iota(p_meta_features_.begin(), p_meta_features_.end(), 0);

        p_labels_ = std::vector<Label>(num_labels_);
        std::iota(p_labels_.begin(), p_labels_.end(), 0);
        if( options_.debug_ ) std::cout << "labels: " << p_labels_.size() << std::endl;

        p_actions_ = std::vector<Action>(num_actions_);
        std::iota(p_actions_.begin(), p_actions_.end(), 0);
        if( options_.debug_ ) std::cout << "actions: " << p_actions_.size() << std::endl;

        p_meta_objects_ = std::vector<MetaObject>(num_meta_objects_);
        std::iota(p_meta_objects_.begin(), p_meta_objects_.end(), 0);
        if( options_.debug_ ) std::cout << "meta-objects: " << p_meta_objects_.size() << std::endl;

        p_unary_ = std::vector<Unary>(num_static_unary_predicates_);
        std::iota(p_unary_.begin(), p_unary_.end(), 0);
        if( options_.debug_ ) std::cout << "unary: " << p_unary_.size() << std::endl;

        p_binary_ = std::vector<Binary>(num_static_binary_predicates_);
        std::iota(p_binary_.begin(), p_binary_.end(), 0);
        if( options_.debug_ ) std::cout << "binary: " << p_binary_.size() << std::endl;


        p_features_ = std::vector<Feature>(num_features_);
        std::iota(p_features_.begin(), p_features_.end(), 0);
        if( options_.debug_ ) std::cout << "features: " << p_features_.size() << std::endl;

        p_states_ = std::vector<State>(num_states_);
        std::iota(p_states_.begin(), p_states_.end(), 0);
        if( options_.debug_ ) std::cout << "states: " << p_states_.size() << std::endl;

        p_transitions_ = std::vector<Transition>(num_transitions_);
        std::iota(p_transitions_.begin(), p_transitions_.end(), 0);
        if( options_.debug_ ) std::cout << "transitions: " << p_transitions_.size() << std::endl;

        p_objects_ = std::vector<Object>(max_num_objects_);
        std::iota(p_objects_.begin(), p_objects_.end(), 0);
        if( options_.debug_ ) std::cout << "objects: " << p_objects_.size() << std::endl;
    }

    // Variables:
    //
    // ==== BEGIN META-LAYER
    //
    // pre0(a,mk) = meta-feature mk appears negated in precondition of action a
    // pre1(a,mk) = meta-feature mk appears positive in precondition of action a
    // eff0(a,mk) = meta-feature mk appears negated in effect of action a
    // eff1(a,mk) = meta-feature mk appears positive in effect of action a
    // label(a,l) = action a is mapped to label l
    //
    // using(mk) = meta-feature mk is used by some action
    // using(a,mk) = meta-feature mk is used by action a
    //
    // arity(p,i) = atom p has arity i
    // atom(mk,p) = meta-feature mk is atom p
    // atom(mk,i,mo) = meta-object mo is mapped to i-th arg of meta-feature mk
    // unary(u,a,mo) = action a uses static unary predicate u on argument mo
    // binary(b,a,mo,mop) = action a uses static binary predicate b on arguments mo and mop
    //
    // args(a,mo) = meta-object mo is *argument* of action a
    // relevant(a,mo,mk,i) = using(a,mk) & atom(mk,i,mo)
    //
    // ==== END META-LAYER
    //
    // map(t,a) = transition t mapped to action a
    // mapf(t,k,mk) = feature k mapped to meta-feature mk for transition t
    // phi(k,s) = value of (boolean) feature k in state s
    //
    // mapt(t,mo,o) = in transition t, meta-object mo is mapped to object o
    // g(k,s,t) = feature k separates states s and t (XOR)
    // free(k,t,a) = feature k isn't affected by transition t (mapped to a)
    //
    // (alt) appl(a,t,s) = action a as in transition t is applicable in state s
    // (alt) mapeq(t,a,tp) = map(t,a) & eq(t,tp)
    // (alt) eq(t,tp) = transitions t and tp are "equivalent"
    // (alt) Z0(t,k,a,s) "=" [ OR { pre0(t,k,a) & mapf(t,k,mk) : mk } => -phi(k,s) ]
    // (alt) Z1(t,k,a,s) "=" [ OR { pre1(t,k,a) & mapf(t,k,mk) : mk } => phi(k,s) ]
    // (alt) X0(a,t,k,mk) "=" [ pre0(a,mk) & mapf(t,k,mk) ]
    // (alt) X1(a,t,k,mk) "=" [ pre1(a,mk) & mapf(t,k,mk) ]
    //
    // ground(k,p) = feature k is ground instance of atom p
    // ground(k,i,o) = i-th arg of feature k is object o
    //
    // W(t,k,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]
    //
    // (enc1) gdiff(k,kp,i) = features k and kp differ in ground atom at i-th arg
    // (enc2) gdiff(k,kp,i) = features k and kp differ in ground atom at i-th arg
    // (enc3) gdiff(k,kp,p) = features k and kp both mapped to atom p
    // (enc3) gdiff(k,kp,p,i) = features k and kp mapped to atom p differ at i-th arg
    //
    // r(l,u,o) = tuple for static unary predicate u in layer l
    // s(l,b,o,op) = tuple for static binary predicate b in layer l
    //
    // U(l,u,a,mo,o) = unary(u,a,mo) & -r(l,u,o)
    // B(l,b,a,mo,mop,o,op) = binary(b,a,mo,mop) & -s(l,b,o,op)
    // ord(o,k,i,s) = ground(k,i,o) & phi(k,s) [ for ordering objects ]
    //
    // FULL (default) ENCODING OF GROUND ACTIONS:
    //
    // gtuple(l,a,<tuple>) = there is (ground) instance of a(<tuple>)
    // G(t,a,<tuple>) = (ground) instance of a(<tuple>) is transition t
    //
    // (alt-default) appl(a,<tuple>,s) = (ground) instance of a(<tuple>) is applicable in state s
    // (alt-default) violated0(a,<tuple>,s,k) = feature k=<some>(<tuple>) is negative precondition of a that doesn't hold in s
    // (alt-default) violated1(a,<tuple>,s,k) = feature k=<some>(<tuple>) is positive precondition of a that doesn't hold in s
    // (alt-default) pre0eq(a,<tuple>,k,mk) => pre0(a,mk) & eq(<tuple>,mk,k)
    // (alt-default) pre1eq(a,<tuple>,k,mk) => pre1(a,mk) & eq(<tuple>,mk,k)
    // (alt-default) eq(<tuple>,mk,k) = meta-feature mk and feature k both mapped to grounded atoms over same <tuple>

    void initialize_variables() override {
        // Meta-layer
        pre0.initialize(*this, "pre0(a,mk)", p_actions_, p_meta_features_);
        pre1.initialize(*this, "pre1(a,mk)", p_actions_, p_meta_features_);
        eff0.initialize(*this, "eff0(a,mk)", p_actions_, p_meta_features_);
        eff1.initialize(*this, "eff1(a,mk)", p_actions_, p_meta_features_);
        label.initialize(*this, "label(a,l)", p_actions_, p_labels_);

        using1.initialize(*this, "using(mk)", p_meta_features_);
        using2.initialize(*this, "using(a,mk)", p_actions_, p_meta_features_);

        arity.initialize(*this, "arity(p,i)", p_atoms_, p_arities_);
        atom1.initialize(*this, "atom(mk,p)", p_meta_features_, p_atoms_);
        atom2.initialize(*this, "atom(mk,i,mo)", p_meta_features_, p_arities_, p_meta_objects_);
        non_static0.initialize(*this, "non_static0(a,mk,p)", p_actions_, p_meta_features_, p_atoms_);
        non_static1.initialize(*this, "non_static1(a,mk,p)", p_actions_, p_meta_features_, p_atoms_);
        unary.initialize(*this, "unary(u,a,mo)", p_unary_, p_actions_, p_meta_objects_);
        binary.initialize(*this, "binary(b,a,mo,mop)", p_binary_, p_actions_, p_meta_objects_, p_meta_objects_);

        args.initialize(*this, "args(a,mo)", p_actions_, p_meta_objects_);
        relevant.initialize(*this, "relevant(a,mo,mk,i)", p_actions_, p_meta_objects_, p_meta_features_, p_arities_);

        // Layers
        map.initialize(*this, "map(t,a)", p_transitions_, p_actions_);
        mapf.initialize(*this, "mapf(t,k,mk)", p_transitions_, p_features_, p_meta_features_);
        phi.initialize(*this, "phi(k,s)", p_features_, p_states_);

        mapt.initialize(*this, "mapt(t,mo,o)", p_transitions_, p_meta_objects_, p_objects_);
        g.initialize(*this, "g(k,s,t)", p_features_, p_states_, p_states_);
        free.initialize(*this, "free(k,t,a)", p_features_, p_transitions_, p_actions_);

        if( encoding("applicable-actions-alt") ) {
            appl.initialize(*this, "appl(a,t,s)", p_actions_, p_transitions_, p_states_);
            mapeq.initialize(*this, "mapeq(t,a,tp)", p_transitions_, p_actions_, p_transitions_);
            eq.initialize(*this, "eq(t,tp)", p_transitions_, p_transitions_);
            Z0.initialize(*this, "Z0(t,k,a,s)", p_transitions_, p_features_, p_actions_, p_states_);
            Z1.initialize(*this, "Z1(t,k,a,s)", p_transitions_, p_features_, p_actions_, p_states_);
            X0.initialize(*this, "X0(a,t,k,mk)", p_actions_, p_transitions_, p_features_, p_meta_features_);
            X1.initialize(*this, "X1(a,t,k,mk)", p_actions_, p_transitions_, p_features_, p_meta_features_);
        }

        ground1.initialize(*this, "ground(k,p)", p_features_, p_atoms_);
        ground2.initialize(*this, "ground(k,i,o)", p_features_, p_arities_, p_objects_);
        W.initialize(*this, "W(t,k,i,mo)", p_transitions_, p_features_, p_arities_, p_meta_objects_);

        if( encoding("features-to-actions-map-alt1") || encoding("features-to-actions-map-alt2") ) {
            gdiff.initialize(*this, "gdiff(k,kp,i)", p_features_, p_features_, p_arities_);
        } else if( encoding("features-to-actions-map-alt3") ) {
            gdiff1.initialize(*this, "gdiff(k,kp,p)", p_features_, p_features_, p_atoms_);
            gdiff2.initialize(*this, "gdiff(k,kp,p,i)", p_features_, p_features_, p_atoms_, p_arities_);
        }

        r.initialize(*this, "r(l,u,o)", p_layers_, p_unary_, p_objects_);
        s.initialize(*this, "s(l,b,o,op)", p_layers_, p_binary_, p_objects_, p_objects_);

        U.initialize(*this, "U(l,u,a,mo,o)", p_layers_, p_unary_, p_actions_, p_meta_objects_, p_objects_);
        B.initialize(*this, "B(l,b,a,mo,mop,o,op)", p_layers_, p_binary_, p_actions_, p_meta_objects_, p_meta_objects_, p_objects_, p_objects_);

        ord.initialize(*this, "ord(o,k,i,s)", p_objects_, p_features_, p_arities_, p_states_);

        // FULL (default) ENCODING OF GROUND ACTIONS:

        // gtuple(l,a,<tuple>) = there is (ground) instance of a(<tuple>)
        gtuple.fill_multipliers(p_layers_, p_actions_);
        for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
            gtuple.fill_multipliers(p_objects_);
        gtuple.initialize_from_multipliers(*this, "gtuple(l,a,<tuple>)");

        // G(t,a,<tuple>) = (ground) instance of a(<tuple>) is transition t
        G.fill_multipliers(p_transitions_, p_actions_);
        for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
            G.fill_multipliers(p_objects_);
        G.initialize_from_multipliers(*this, "G(t,a,<tuple>)");

        if( encoding("applicable-actions-tuples") ) {
            // appl(a,<tuple>,s) = (ground) instance of a(<tuple>) is applicable in state s
            appl2.fill_multipliers(p_actions_);
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                appl2.fill_multipliers(p_objects_);
            appl2.fill_multipliers(p_states_);
            appl2.initialize_from_multipliers(*this, "appl(a,<tuple>,s)");

            // violated0(a,<tuple>,s,k) = feature k=<some>(<tuple>) is negative precondition of a that doesn't hold in s
            violated0.fill_multipliers(p_actions_);
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                violated0.fill_multipliers(p_objects_);
            violated0.fill_multipliers(p_states_);
            violated0.fill_multipliers(p_features_);
            violated0.initialize_from_multipliers(*this, "violated0(a,<tuple>,s,k)");

            // violated1(a,<tuple>,s,k) = feature k=<some>(<tuple>) is positive precondition of a that doesn't hold in s
            violated1.fill_multipliers(p_actions_);
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                violated1.fill_multipliers(p_objects_);
            violated1.fill_multipliers(p_states_);
            violated1.fill_multipliers(p_features_);
            violated1.initialize_from_multipliers(*this, "violated1(a,<tuple>,s,k)");

            // pre0eq(a,<tuple>,k,mk) => pre0(a,mk) & eq(<tuple>,mk,k)
            pre0eq.fill_multipliers(p_actions_);
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                pre0eq.fill_multipliers(p_objects_);
            pre0eq.fill_multipliers(p_features_);
            pre0eq.fill_multipliers(p_meta_features_);
            pre0eq.initialize_from_multipliers(*this, "pre0eq(a,<tuple>,k,mk)");

            // pre1eq(a,<tuple>,k,mk) => pre1(a,mk) & eq(<tuple>,mk,k)
            pre1eq.fill_multipliers(p_actions_);
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                pre1eq.fill_multipliers(p_objects_);
            pre1eq.fill_multipliers(p_features_);
            pre1eq.fill_multipliers(p_meta_features_);
            pre1eq.initialize_from_multipliers(*this, "pre1eq(a,<tuple>,k,mk)");

            // eq(<tuple>,mk,k) = meta-feature mk and feature k both mapped to grounded atoms over same <tuple>
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                eq2.fill_multipliers(p_objects_);
            eq2.fill_multipliers(p_meta_features_);
            eq2.fill_multipliers(p_features_);
            eq2.initialize_from_multipliers(*this, "eq(<tuple>,mk,k)");
        }
    }

    // formula groups
    void start_group(const std::string &group_desc) {
        group_descriptions_.emplace_back(group_desc);
        imp_offsets_.emplace_back(num_implications(), group_desc);
        add_comment(group_desc);
    }
    void end_group(std::ostream &os) {
        os << "(" << group_descriptions_.size() << ") " << group_descriptions_.back()
           << ": #implications="
           << num_implications() - imp_offsets_.back().first
           << std::endl;
    }
    void end_group(std::ostream &os, const std::string &alt_desc, int n) {
        os << "(" << group_descriptions_.size() << ") " << alt_desc
           << ": #implications=" << n << std::endl;
    }


    // Base formulas:
    //
    // ==== BEGIN META-LAYER
    //
    // Basic definitions:
    //  (1) using(mk) <=> OR { using(a,mk) : a }
    //  (2) using(a,mk) <=> pre0(a,mk) v pre1(a,mk) v eff0(a,mk) v eff1(a,mk)
    //
    // Consistent preconditions, effects, and labeling:
    //  (3) -pre0(a,mk) v -pre1(a,mk)
    //  (4) -eff0(a,mk) v -eff1(a,mk)
    //  (5) AT-MOST-1 { label(a,l) : l }
    //
    // Effects are non-redundant:
    //  (6) eff0(a,mk) => -pre0(a,mk)
    //  (7) eff1(a,mk) => -pre1(a,mk)
    //
    // Disjoint meta-features in actions:
    //  (8) (optional) AT-MOST-1 { using(a,mk) : a }
    //
    // Consistent arities for atoms:
    //  (9) EXACTLY-1 { arity(p,i) : 0 <= i <= max-arity }
    //
    // Mapping of meta-features to atoms:
    // (10) EXACTLY-1 { atom(mk,p) : p }
    // (11) AT-MOST-1 { atom(mk,i,mo) : mo }
    // (12) atom(mk,p) & atom(mk,i,mo) => OR { arity(p,j) : i <= j <= max-arity }
    // (13) atom(mk,p) & arity(p,i) => OR { atom(mk,j,mo) : mo } [ 1 <= j <= i ]
    // (14) atom(mk,p) & arity(p,i) => -atom(mk,j,mo) [ i < j ]
    // (15) -atom(mk,0,mo)
    //
    // 1-1 mapping of meta-features to atoms:
    // (16) <strict-lexicographic-ordering-meta-features>
    //
    // Atoms must be non-static:
    // (xx) OR_{a,mk} non-static0(a,mk,p) v non-static1(a,mk,p)
    // (xx) non-static0(a,mk,p) => atom(mk,p) & pre1(a,mk) & eff0(a,mk)
    // (xx) non-static1(a,mk,p) => atom(mk,p) & pre0(a,mk) & eff1(a,mk)
    //
    // Relevant arguments for actions:
    // (17) using(a,mk) & atom(mk,i,mo) => args(a,mo)
    // (18) args(a,mo) => OR { relevant(a,mo,mk,i) : mk, i }
    // (19) relevant(a,mo,mk,i) <=> using(a,mk) & atom(mk,i,mo)
    //
    // Arity for actions and atoms:
    // (20) -args(a,mo_i) [ i >= <arity-action-a> ]
    // (21) -arity(p,i) [ i > <arity-atom-p> ]
    //
    // (Regularizer) Exact arities:
    // (20x) args(a,mo_i) [ 0 <= i < <arity-action-a> ]
    // (21x) arity(p,i) [ i == <arity-atom-p> ]
    //
    // Static predicates on relevant arguments:
    // (22) unary(u,a,mo) => args(a,mo)
    // (23) binary(b,a,mo,mop) => args(a,mo) & args(a,mop)
    //
    // ==== END META-LAYER
    //
    // Consistent mapping of transitions to actions:
    // (24) EXACTLY-1 { map(t,a) : a }
    //
    // Consistent mapping of features to meta-features in actions:
    // (25) AT-MOST-1 { mapf(t,k,mk) : mk }
    // (26) AT-MOST-1 { mapf(t,k,mk) : k }
    //
    // Consistency between map, mapf, labeling, and using:
    // (27) map(t,a) => label(a,t.label)
    // (28) map(t,a) & mapf(t,k,mk) => using(a,mk)
    // (29) map(t,a) & using(a,mk) => OR { mapf(t,k,mk) : k }
    //
    // Definition of free(k,t,a):
    // (30) map(t,a) & AND { -mapf(t,k,mk) : mk } => free(k,t,a)
    // (31) map(t,a) & mapf(t,k,mk) => [ -eff0(a,mk) & -eff1(a,mk) <=> free(k,t,a) ]
    //
    // Transitions:
    // (32) <transitions>
    // (32.1) map(t,a) & mapf(t,k,mk) & pre0(a,mk) => -phi(k,t.src)
    // (32.2) map(t,a) & mapf(t,k,mk) & pre1(a,mk) => phi(k,t.src)
    // (32.3) map(t,a) & mapf(t,k,mk) & eff0(a,mk) => -phi(k,t.dst)
    // (32.4) map(t,a) & mapf(t,k,mk) & eff1(a,mk) => phi(k,t.dst)
    //
    // Inertia:
    // (33) map(t,a) => [ free(k,t,a) <=> [ phi(k,t.src) <=> phi(k,t.dst) ] ]
    //
    // Applicable actions must be applied:
    // (34) (alt) appl(a,t,s) => OR { mapeq(tp,a,t) : tp.src = s, tp.label = t.label } [ t.src != s ]
    // (35) (alt) mapeq(t,a,tp) => map(t,a) & eq(t,tp)
    // (36) (alt) eq(t,tp) => [ mapf(t,k,mk) <=> mapf(tp,k,mk) ]
    //
    // Definition of g(k,s,t):
    // (37) <def g(k,s,t)>
    // (37.1) g(k,s,t) => phi(k,s) v phi(k,t)
    // (37.2) g(k,s,t) => -phi(k,s) v -phi(k,t)
    // (37.3) phi(k,s) & -phi(k,t) => g(k,s,t)
    // (37.4) -phi(k,s) & phi(k,t) => g(k,s,t)
    //
    // Separate different states using features:
    // (38) OR { g(k,s,t) : k } [ s < t ]
    //
    // Definition of appl(a,t,s):
    // (39) (alt) map(t,a) & AND { Z0(t,k,a,s) : k } & AND { Z1(t,k,a,s) : k } => appl(a,t,s)
    //
    // Definition of Z0(t,k,a,s):
    // (40) (alt) -phi(k,s) => Z0(t,k,a,s)
    // (41) (alt) AND { -X0(a,t,k,mk) : mk } => Z0(t,k,a,s)
    // (42) (alt) X0(a,t,k,mk) => pre0(a,mk) & mapf(t,k,mk)
    //
    // Definition of Z1(t,k,a,s):
    // (43) (alt) phi(k,s) => Z1(t,k,a,s)
    // (44) (alt) AND { -X1(a,t,k,mk) : mk } => Z1(t,k,a,s)
    // (45) (alt) X1(a,t,k,mk) => pre1(a,mk) & mapf(t,k,mk)
    //
    // Mapping of features to grounded atoms:
    // (46) EXACTLY-1 { ground(k,p) : p }
    // (47) AT-MOST-1 { ground(k,i,o) : o }
    // (48) ground(k,p) & ground(k,i,o) => OR { arity(p,j) : i <= j <= max-arity }
    // (49) ground(k,p) & arity(p,i) => OR { ground(k,j,o) : o } [ 1 <= j <= i ]
    // (50) ground(k,p) & arity(p,i) => -ground(k,j,o) [ i < j ]
    // (51) -ground(k,0,o)
    //
    // 1-1 mapping of features to grounded atoms:
    // (52.1.1) (enc1) ground(k,p) & ground(kp,p) & arity(p,i) => OR { gdiff(k,kp,j) : 1 <= j <= i }
    // (52.1.2) (enc1) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
    // (52.1.3) (enc1) gdiff(k,kp,i) => OR { ground(k,i,o) : o } & OR { ground(kp,i,o) : o }
    // (52.1.4) (enc1) gdiff(k,kp,i) => AND { -gdiff(k,kp,j) : 1 <= j < i }
    //
    // (52.2.1) (enc2) arity(p,0) => -ground(k,p) v -ground(kp,p)
    // (52.2.2) (enc2) OR { gdiff(k,kp,i) : 0 <= i <= max-arity }
    // (52.2.3) (enc2) gdiff(k,kp,0) => -ground(k,p) v -ground(kp,p)
    // (52.2.4) (enc2) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
    // (52.2.5) (enc2) gdiff(k,kp,i) => OR { ground(k,i,o) : o } v OR { ground(kp,i,o) : o }
    // (52.2.6) (enc2) [DISABLED] AT-MOST-1 { gdiff(k,kp,i) : 0 <= i <= max-arity }
    //
    // (52.3.1) (enc3) gdiff(k,kp,p) <=> ground(k,p) & ground(kp,p)
    // (52.3.2) (enc3) gdiff(k,kp,p) & arity(p,i) => OR { gdiff(k,kp,p,j) : 1 <= j <= i }
    // (52.3.3) (enc3) gdiff(k,kp,p,j) => gdiff(k,kp,p) & OR { arity(p,i) : j <= i <= max-arity }
    // (52.3.4) (enc3) gdiff(k,kp,p,i) => -ground(k,i,o) v -ground(kp,i,o)
    //
    // (52.4.1) (enc4) implemented by <strict-lexicographic-ordering-features-in-layers>
    //
    // 1-1 mapping of features to grounded atoms:
    // (53.1) (enc123) <non-strict-lexicographic-ordering-features-in-layers>
    // (53.2) (enc4) <strict-lexicographic-ordering-features-in-layers>
    //
    // Consistent mapping between features and meta-features
    // (54) mapf(t,k,mk) => [ atom(mk,p) <=> ground(k,p) ]
    // (55) mapf(t,k,mk) & atom(mk,i,mo) => OR { ground(k,i,o) : o }
    // (56) mapf(t,k,mk) & ground(k,i,o) => OR { atom(mk,i,mo) : mo }
    //
    // Definition of U(l,u,a,mo,o) and B(l,b,a,mo,mop,o,op):
    // (57) U(l,u,a,mo,o) <=> unary(u,a,mo) & -r(l,u,o)
    // (58) B(l,b,a,mo,mop,o,op) <=> binary(b,a,mo,mop) & -s(l,b,o,op)
    //
    // Definition of mapt(t,mo,o):
    // (59) AT-MOST-1 { mapt(t,mo,o) : o }
    // (60) map(t,a) & args(a,mo) => OR { mapt(t,mo,o) : o }
    // (61) map(t,a) & mapt(t,mo,o) => args(a,mo)
    //
    // Cross consistency between schemas and transitions:
    // (62) mapf(t,k,mk) & atom(mk,i,mo) => W(t,k,i,mo)
    // (63) W(t,k,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]
    //
    // One-hot encoding of variables:
    // (64) <bounds-feature-sum-layers>
    // (65) (enc123) -phi(k,s0) [ k < #features - <upper-bound> ]
    // (66) (enc123) phi(k,s0) [ k >= #features - <lower-bound> ]
    //
    // Explanation of existing grounded actions (common to different encodings):
    // (67) map(t,a) & mapt(t,mo,o) & unary(u,a,mo) => r(t.layer,u,o)
    // (68) map(t,a) & mapt(t,mo,o) & mapt(t,mop,op) & binary(b,a,mo,mop) => s(t.layer,b,o,op)
    //
    // (Symmetries) Ordered arities for atoms:
    // (69) arity(p,i) => OR { arity(p-1,j) : 0 <= j <= i }
    //
    // (Symmetries) Ordered use of action arguments:
    // (70) args(a,mo) => args(a,mop) [ 1 + mop = mo ]
    //
    // (Symmetries) Definition of ord(o,k,i,s):
    // (71) ord(o,k,i,s) <=> ground(k,i,o) & phi(k,s)
    //
    // (Symmetries) Ordered objects at each layer
    // (72) <non-strict-lexicographic-ordering-objects-in-layers>
    //
    // FULL (default) ENCODING OF GROUND ACTIONS:
    //
    // Explanation of non-existing grounded actions:
    // (73) -gtuple(l,a,<tuple>) => OR { -args(a,mo_i) : o_i > 0 } v
    //                              OR { U(l,u,a,mo_i,o_i) : u, i } v
    //                              OR { B(l,b,a,mo_i,mo_j,o_i,o_j) : b, i < j }
    //
    // Definition of gtuple(l,a,<tuple>) and G(t,a,<tuple>):
    // (74) gtuple(l,a,<tuple>) => OR { G(t,a,<tuple>) : t.layer = l }
    // (75) G(t,a,<tuple>) => gtuple(t.layer,a,<tuple>)
    // (76) G(t,a,<tuple>) => map(t,a) &
    //                        AND { mapt(t,mo_i,o_i) : o_i > 0 }
    //                        AND { args(a,mo_i) => mapt(t,mo_i,o_i) : o_i = 0 }
    //
    // Explanation of existing grounded actions:
    // (77) AT-MOST-1 { G(t,a,<tuple>) : t.src = s } [ s, a(<tuple>) ]
    // (78) EXACTLY-1 { G(t,a,<tuple>) : a, <tuple> } [ t ]
    // (79) [DISABLED] AT-LEAST-1 { G(t,a,<tuple>) : a, <tuple> } [ t ]
    // (80) [SUBSUMED BY 62] gtuple(l,a,<tuple>) & unary(u,a,moi) => r(l,u,oi)
    // (81) [SUBSUMED BY 63] gtuple(l,a,<tuple>) & binary(b,a,moi,moj) => s(l,b,oi,oj)
    //
    // Applicable actions must be applied:
    // (82) (tuples-default) G(t,a,<tuple>) => appl(a,<tuple>,t.src)
    // (83) (tuples-default) appl(a,<tuple>,s) => OR { G(t,a,<tuple>) : t.src = s }
    // (84) (tuples-default) -appl(a,<tuple>,s) => -gtuple(s.layer,a,<tuple>) v
    //                                             OR { violated0(a,<tuple>,s,k) : k } v
    //                                             OR { violated1(a,<tuple>,s,k) : k }
    // (85) (tuples-default) violated0(a,<tuple>,s,k) => phi(k,s) & OR { pre0eq(a,<tuple>,k,mk) : mk }
    // (86) (tuples-default) violated1(a,<tuple>,s,k) => -phi(k,s) & OR { pre1eq(a,<tuple>,k,mk) : mk }
    // (87) (tuples-default) pre0eq(a,<tuple>,k,mk) => pre0(a,mk) & eq(<tuple>,mk,k)
    // (88) (tuples-default) pre1eq(a,<tuple>,k,mk) => pre1(a,mk) & eq(<tuple>,mk,k)
    // (89) (tuples-default) eq(<tuple>,mk,k) => [ atom(mk,p) <=> ground(k,p) ]
    // (90) (tuples-default) eq(<tuple>,mk,k) & atom(mk,i,mo_j) => ground(k,i,o_j)

    // Basic definitions:
    // using(mk) <=> OR { using(a,mk) : a }
    void build_formulas_using1_iff_using2(std::ostream &os) {
        start_group("using(mk) <=> OR { using(a,mk) : a }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            assert((0 <= mk) && (mk < num_meta_features_));

            // using(mk) => OR { using(a,mk) : a }
            SAT::Implication IP({ 1 + using1(mk) }, { });
            for( Action a = Action(0); a < num_actions_; ++a )
                IP.add_consequent(1 + using2(a, mk));
            add_implication(IP);

            // using(mk) <= OR { using(a,mk) : a }
            for( Action a = Action(0); a < num_actions_; ++a )
                add_implication({ 1 + using2(a, mk) }, { 1 + using1(mk) });
        };

        using1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Basic definitions:
    // using(a,mk) <=> pre0(a,mk) v pre1(a,mk) v eff0(a,mk) v eff1(a,mk)
    void build_formulas_using2_iff_pre0_pre1_eff0_eff1(std::ostream &os) {
        start_group("using(a,mk) <=> pre0(a,mk) v pre1(a,mk) v eff0(a,mk) v eff1(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            // using(a,mk) => pre0(a,mk) v pre1(a,mk) v eff0(a,mk) v eff1(a,mk)
            // pre0(a,mk) => using(a,mk)
            // pre1(a,mk) => using(a,mk)
            // eff0(a,mk) => using(a,mk)
            // eff1(a,mk) => using(a,mk)
            add_implication({ 1 + using2(a, mk) }, { 1 + pre0(a, mk), 1 + pre1(a, mk), 1 + eff0(a, mk), 1 + eff1(a, mk) });
            add_implication({ 1 + pre0(a, mk) }, { 1 + using2(a, mk) });
            add_implication({ 1 + pre1(a, mk) }, { 1 + using2(a, mk) });
            add_implication({ 1 + eff0(a, mk) }, { 1 + using2(a, mk) });
            add_implication({ 1 + eff1(a, mk) }, { 1 + using2(a, mk) });
        };

        using2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent preconditions, effects, and labeling:
    // -pre0(a,mk) v -pre1(a,mk)
    void build_formulas_then_pre0_pre1(std::ostream &os) {
        start_group("-pre0(a,mk) v -pre1(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            // -pre0(a,mk) v -pre1(a,mk)
            add_implication({ }, { -(1 + pre0(a, mk)), -(1 + pre1(a, mk)) });
        };

        pre0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent preconditions, effects, and labeling:
    // -eff0(a,mk) v -eff1(a,mk)
    void build_formulas_then_eff0_eff1(std::ostream &os) {
        start_group("-eff0(a,mk) v -eff1(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            // -eff0(a,mk) v -eff1(a,mk)
            add_implication({ }, { -(1 + eff0(a, mk)), -(1 + eff1(a, mk)) });
        };

        eff0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent preconditions, effects, and labeling:
    // AT-MOST-1 { label(a,l) : l }
    void build_formulas_at_most_1_label(std::ostream &os) {
        start_group("AT-MOST-1 { label(a,l) : l }");

        for( Action a = Action(0); a < num_actions_; ++a ) {
            std::vector<int> literals(num_labels_);
            for( Label l = Label(0); l < num_labels_; ++l )
                literals[l] = 1 + label(a, l);
            at_most_1(std::string("label-at-most-1(") + std::to_string(a) + ")", literals);
        }

        end_group(os);
    }

    // Effects are non-redundant:
    // eff0(a,mk) => -pre0(a,mk)
    // eff1(a,mk) => -pre1(a,mk)
    void build_formulas_eff_then_pre(std::ostream &os) {
        start_group("eff0(a,mk) => -pre0(a,mk) and eff1(a,mk) => -pre1(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            // eff0(a,mk) => -pre0(a,mk)
            // eff1(a,mk) => -pre1(a,mk)
            add_implication({ 1 + eff0(a, mk) }, { -(1 + pre0(a, mk)) });
            add_implication({ 1 + eff1(a, mk) }, { -(1 + pre1(a, mk)) });
        };

        eff0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (Regularizer) Disjoint meta-features in actions:
    // AT-MOST-1 { using(a,mk) : a }
    void build_formulas_at_most_1_using2(std::ostream &os) {
        start_group("AT-MOST-1 { using(a,mk) : a }");

        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            std::vector<int> literals(num_actions_);
            for( Action a = Action(0); a < num_actions_; ++a )
                literals[a] = 1 + using2(a, mk);
            at_most_1(std::string("using2-at-most-1(") + std::to_string(mk) + ")", literals);
        }

        end_group(os);
    }

    // Consistent arities for atoms:
    // EXACTLY-1 { arity(p,i) : 0 <= i <= max-arity }
    void build_formulas_exactly_1_arity(std::ostream &os) {
        start_group("EXACTLY-1 { arity(p,i) : 0 <= i <= max-arity }");

        for( Atom p = Atom(0); p < num_atoms_; ++p ) {
            std::vector<int> literals(1 + max_arity_);
            for( Arity i = Arity(0); i < 1 + max_arity_; ++i )
                literals[i] = 1 + arity(p, i);
            exactly_1(std::string("arity-exactly-1(") + std::to_string(p) + ")", literals);
        }

        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // EXACTLY-1 { atom(mk,p) : p }
    void build_formulas_exactly_1_atom1(std::ostream &os) {
        start_group("EXACTLY-1 { atom(mk,p) : p }");

        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            std::vector<int> literals(num_atoms_);
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                literals[p] = 1 + atom1(mk, p);
            exactly_1(std::string("atom-exactly-1(") + std::to_string(mk) + ")", literals);
        }

        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // AT-MOST-1 { atom(mk,i,mo) : mo }
    void build_formulas_at_most_1_atom2(std::ostream &os) {
        start_group("AT-MOST-1 { atom(mk,i,mo) : mo }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            Arity i(tuple.at(1));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( i > 0 ) {
                std::vector<int> literals(num_meta_objects_);
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    literals[mo] = 1 + atom2(mk, i, mo);
                at_most_1(std::string("atom-at-most-1(") + std::to_string(mk) + "," + std::to_string(i) + ")", literals);
            }
        };

        SAT::VarSet tmp(p_meta_features_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // atom(mk,p) & atom(mk,i,mo) => OR { arity(p,j) : i <= j <= max-arity }
    void build_formulas_atom1_atom2_then_arity(std::ostream &os) {
        start_group("atom(mk,p) & atom(mk,i,mo) => OR { arity(p,j) : i <= j }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2));
            MetaObject mo(tuple.at(3));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( i > 0 ) {
                // atom(mk,p) & atom(mk,i,mo) => OR { arity(p,j) : i <= j <= max-arity }
                SAT::Implication IP({ 1 + atom1(mk, p), 1 + atom2(mk, i, mo) }, { });
                for( Arity j = i; j < 1 + max_arity_; ++j )
                    IP.add_consequent(1 + arity(p, j));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_meta_features_, p_atoms_, p_arities_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // atom(mk,p) & arity(p,i) => OR { atom(mk,j,mo) : mo } [ 1 <= j <= i ]
    void build_formulas_atom1_arity_then_atom2(std::ostream &os) {
        start_group("atom(mk,p) & arity(p,i) => OR { atom(mk,j,mo) : mo }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2)), j(tuple.at(3));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= j) && (j < 1 + max_arity_));

            if( (0 < j) && (j <= i) ) {
                // atom(mk,p) & arity(p,i) => OR { atom(mk,j,mo) : mo } [ 1 <= j <= i ]
                SAT::Implication IP({ 1 + atom1(mk, p), 1 + arity(p, i) }, { });
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    IP.add_consequent(1 + atom2(mk, j, mo));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_meta_features_, p_atoms_, p_arities_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // atom(mk,p) & arity(p,i) => -atom(mk,j,mo) [ i < j ]
    void build_formulas_atom1_arity_then_atom2_2(std::ostream &os) {
        start_group("atom(mk,p) & arity(p,i) => -atom(mk,j,mo) [i < j]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2)), j(tuple.at(3));
            MetaObject mo(tuple.at(4));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= j) && (j < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( i < j ) {
                // atom(mk,p) & arity(p,i) => -atom(mk,j,mo) [ i < j ]
                add_implication({ 1 + atom1(mk, p), 1 + arity(p, i) }, { -(1 + atom2(mk, j, mo)) });
            }
        };

        SAT::VarSet tmp(p_meta_features_, p_atoms_, p_arities_, p_arities_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of meta-features to atoms:
    // -atom(mk,0,mo)
    void build_formulas_atom2(std::ostream &os) {
        start_group("-atom(mk,0,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            MetaFeature mk(tuple.at(0));
            MetaObject mo(tuple.at(1));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            // -atom(mk,0,mo)
            add_unit(-(1 + atom2(mk, 0, mo)));
        };

        SAT::VarSet tmp(p_meta_features_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of meta-features to atoms:
    // <strict-lexicographic-ordering-meta-features>
    void build_formulas_lex_ordering_meta_layer(std::ostream &os) {
        start_group("<strict-lexicographic-ordering-meta-features>");

        std::string prefix("meta-layer");
        std::vector<std::vector<int> > lvectors;

        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            std::vector<int> meta_feature;

            // used meta-features appear first
            meta_feature.push_back(-(1 + using1(mk)));

            if( regularizer("disjoint-meta-features") ) {
                // order by usage in actions
                for( Action a = Action(0); a < num_actions_; ++a )
                    meta_feature.push_back(-(1 + using2(a, mk)));
            }

            // order by atom names
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                meta_feature.push_back(-(1 + atom1(mk, p)));

            // order by arguments
            for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    meta_feature.push_back(-(1 + atom2(mk, i, mo)));
            }

            // wrong as it may lead to different mks mapped to same action when the mks appear in different actions
            if( false ) { // && !options_.disjoint_meta_features_ ) {
                // order by usage in actions
                for( Action a = Action(0); a < num_actions_; ++a )
                    meta_feature.push_back(-(1 + using2(a, mk)));
            }

            lvectors.emplace_back(std::move(meta_feature));
        }

        // lex ordering in meta-layer must be weak!
        if( !lvectors.empty() ) lex_ordering(prefix, lvectors, !options_.weak_lex_orderings_);

        end_group(os);
    }

    // Atoms must be non-static:
    // OR { non-static0(a,mk,p) v non_static1(a,mk,p) : a, mk }
    void build_formulas_non_static0_non_static1(std::ostream &os) {
        start_group("OR { non-static0(a,mk,p) v non_static1(a,mk,p) : a, mk }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Atom p(tuple.at(0));
            assert((0 <= p) && (p < num_atoms_));

            // OR { non-static0(a,mk,p) v non_static1(a,mk,p) : a, mk }
            SAT::Implication IP;
            for( Action a = Action(0); a < num_actions_; ++a ) {
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
                    IP.add_consequent(1 + non_static0(a, mk, p));
                    IP.add_consequent(1 + non_static1(a, mk, p));
                }
            }
            add_implication(IP);
        };

        SAT::VarSet tmp(p_atoms_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Atoms must be non-static:
    // non-static0(a,mk,p) => atom(mk,p) & pre1(a,mk) & eff0(a,mk)
    void build_formulas_non_static0_then_atom_pre1_eff0(std::ostream &os) {
        start_group("non-static0(a,mk,p) => atom(mk,p) & pre1(a,mk) & eff0(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            Atom p(tuple.at(2));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= p) && (p < num_atoms_));

            // non-static0(a,mk,p) => atom(mk,p)
            // non-static0(a,mk,p) => pre1(a,mk)
            // non-static0(a,mk,p) => eff0(a,mk)
            add_implication({ 1 + non_static0(a, mk, p) }, { 1 + atom1(mk, p) });
            add_implication({ 1 + non_static0(a, mk, p) }, { 1 + pre1(a, mk) });
            add_implication({ 1 + non_static0(a, mk, p) }, { 1 + eff0(a, mk) });
        };

        non_static0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Atoms must be non-static:
    // non-static1(a,mk,p) => atom(mk,p) & pre0(a,mk) & eff1(a,mk)
    void build_formulas_non_static1_then_atom_pre0_eff1(std::ostream &os) {
        start_group("non-static1(a,mk,p) => atom(mk,p) & pre0(a,mk) & eff1(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            Atom p(tuple.at(2));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= p) && (p < num_atoms_));

            // non-static1(a,mk,p) => atom(mk,p)
            // non-static1(a,mk,p) => pre0(a,mk)
            // non-static1(a,mk,p) => eff1(a,mk)
            add_implication({ 1 + non_static1(a, mk, p) }, { 1 + atom1(mk, p) });
            add_implication({ 1 + non_static1(a, mk, p) }, { 1 + pre0(a, mk) });
            add_implication({ 1 + non_static1(a, mk, p) }, { 1 + eff1(a, mk) });
        };

        non_static1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Relevant arguments for actions:
    // using(a,mk) & atom(mk,i,mo) => args(a,mo)
    void build_formulas_using2_atom2_then_args(std::ostream &os) {
        start_group("using(a,mk) & atom(mk,i,mo) => args(a,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( 0 < i ) {
                // using(a,mk) & atom(mk,i,mo) => args(a,mo)
                add_implication({ 1 + using2(a, mk), 1 + atom2(mk, i, mo) }, { 1 + args(a, mo) });
            }
        };

        SAT::VarSet tmp(p_actions_, p_meta_objects_, p_meta_features_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Relevant arguments for actions:
    // args(a,mo) => OR { relevant(a,mo,mk,i) : mk, i }
    void build_formulas_args_then_relevant(std::ostream &os) {
        start_group("args(a,mo) => OR { relevant(a,mo,mk,i) : mk, i }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            // args(a,mo) => OR { relevant(a,mo,mk,i) : mk, i }
            SAT::Implication IP({ 1 + args(a, mo) }, { });
            for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
                for( Arity i = Arity(1); i < 1 + max_arity_; ++i )
                    IP.add_consequent(1 + relevant(a, mo, mk, i));
            }
            add_implication(IP);
        };

        SAT::VarSet tmp(p_actions_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Relevant arguments for actions:
    // relevant(a,mo,mk,i) <=> using(a,mk) & atom(mk,i,mo)
    void build_formulas_relevant_iff_using2_atom2(std::ostream &os) {
        start_group("relevant(a,mo,mk,i) <=> using(a,mk) & atom(mk,i,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( 0 < i ) {
                // relevant(a,mo,mk,i) => using(a,mk)
                // relevant(a,mo,mk,i) => atom(mk,i,mo)
                // relevant(a,mo,mk,i) <= using(a,mk) & atom(mk,i,mo)
                add_implication({ 1 + relevant(a, mo, mk, i) }, { 1 + using2(a, mk) });
                add_implication({ 1 + relevant(a, mo, mk, i) }, { 1 + atom2(mk, i, mo) });
                add_implication({ 1 + using2(a, mk), 1 + atom2(mk, i, mo) }, { 1 + relevant(a, mo, mk, i) });
            }
        };

        relevant.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Arity for actions and atoms:
    // -args(a,mo_i) [ i >= <arity-action-a> ]
    void build_formulas_args(std::ostream &os) {
        start_group("-args(a,mo_i) [ i >= <arity-action-a> ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            //assert(!regularizer("disjoint-meta-features"));

            if( arity_for_actions_.at(a) <= mo )
                add_unit(-(1 + args(a, mo)));
        };

        args.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Arity for actions and atoms:
    // -arity(p,i) [ i > <arity-atom-p> ]
    void build_formulas_arity(std::ostream &os) {
        start_group("-arity(p,i) [ i > <arity-atom-p> ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Atom p(tuple.at(0));
            Arity i(tuple.at(1));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));
            if( arity_for_atoms_.at(p) < i )
                add_unit(-(1 + arity(p, i)));
        };

        arity.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Arity for actions and atoms:
    // args(a,mo_i) [ 0 <= i < <arity-action-a> ]
    void build_formulas_exact_args(std::ostream &os) {
        start_group("args(a,mo_i) [ 0 <= i < <arity-action-a> ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            //assert(!regularizer("disjoint-meta-features"));

            if( mo < arity_for_actions_.at(a) )
                add_unit(1 + args(a, mo));
        };

        args.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Arity for actions and atoms:
    // arity(p,i) [ i == <arity-atom-p> ]
    void build_formulas_exact_arity(std::ostream &os) {
        start_group("arity(p,i) [ i == <arity-atom-p> ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Atom p(tuple.at(0));
            Arity i(tuple.at(1));
            assert((0 <= p) && (p < num_atoms_));
            if( i == arity_for_atoms_.at(p) )
                add_unit(1 + arity(p, i));
        };

        arity.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Static predicates on relevant arguments:
    // unary(u,a,mo) => args(a,mo)
    void build_formulas_unary_then_args(std::ostream &os) {
        start_group("unary(u,a,mo) => args(a,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Unary u(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2));
            assert((0 <= u) && (u < num_static_unary_predicates_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            // unary(u,a,mo) => args(a,mo)
            add_implication({ 1 + unary(u, a, mo) }, { 1 + args(a, mo) });
        };

        unary.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Static predicates on relevant arguments:
    // binary(b,a,mo,mop) => args(a,mo) & args(a,mop)
    void build_formulas_binary_then_args_args(std::ostream &os) {
        start_group("binary(b,a,mo,mop) => args(a,mo) & args(a,mop)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Binary b(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2)), mop(tuple.at(3));
            assert((0 <= b) && (b < num_static_binary_predicates_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= mop) && (mop < num_meta_objects_));

            // binary(b,a,mo,mop) => args(a,mo)
            // binary(b,a,mo,mop) => args(a,mop)
            add_implication({ 1 + binary(b, a, mo, mop) }, { 1 + args(a, mo) });
            add_implication({ 1 + binary(b, a, mo, mop) }, { 1 + args(a, mop) });
        };

        binary.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent mapping of transitions to actions:
    // EXACTLY-1 { map(t,a) : a }
    void build_formulas_exactly_1_map(std::ostream &os) {
        start_group("EXACTLY-1 { map(t,a) : a }");

        for( Transition t = Transition(0); t < num_transitions_; ++t ) {
            std::vector<int> literals(num_actions_);
            for( Action a = Action(0); a < num_actions_; ++a )
                literals[a] = 1 + map(t, a);
            exactly_1(std::string("map-exactly-1(") + std::to_string(t) + ")", literals);
        }

        end_group(os);
    }

    // Consistent mapping of features to meta-features in actions:
    // AT-MOST-1 { mapf(t,k,mk) : mk }
    void build_formulas_at_most_1_mapf_by_k(std::ostream &os) {
        start_group("AT-MOST-1 { mapf(t,k,mk) : mk }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

            if( layer == f_layer_.at(k) ) {
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                std::vector<int> literals(num_meta_features_);
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk )
                    literals[mk] = 1 + mapf(t, k, mk);
                at_most_1(std::string("mapf-k-at-most-1(") + std::to_string(t) + "," + std::to_string(k) + ")", literals);
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent mapping of features to meta-features in actions:
    // AT-MOST-1 { mapf(t,k,mk) : k }
    void build_formulas_at_most_1_mapf_by_mk(std::ostream &os) {
        start_group("AT-MOST-1 { mapf(t,k,mk) : k }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= mk) && (mk < num_meta_features_));

            std::vector<int> literals(features_per_layer_[layer].size());
            for( int i = 0; i < int(features_per_layer_[layer].size()); ++i )
                literals[i] = 1 + mapf(t, features_per_layer_[layer].at(i), mk);
            at_most_1(std::string("mapf-mk-at-most-1(") + std::to_string(t) + "," + std::to_string(mk) + ")", literals);
        };

        SAT::VarSet tmp(p_transitions_, p_meta_features_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistency between map, mapf, labeling, and using:
    // map(t,a) => label(a,t.label)
    void build_formulas_map_then_label(std::ostream &os) {
        start_group("map(t,a) => label(a,t.label)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Action a(tuple.at(1));
            assert((0 <= t) && (t < num_transitions_));
            assert((0 <= a) && (a < num_actions_));

            // map(t,a) => label(a,l)
            Label l = tr_label_.at(t);
            add_implication({ 1 + map(t, a) }, { 1 + label(a, l) });
        };

        map.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistency between map, mapf, labeling, and using:
    // map(t,a) & mapf(t,k,mk) => using(a,mk)
    void build_formulas_map_mapf_then_using2(std::ostream &os) {
        start_group("map(t,a) & mapf(t,k,mk) => using(a,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Action a(tuple.at(3));
            Layer layer = tr_layer_.at(t);
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

            if( layer == f_layer_.at(k) ) {
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // map(t,a) & mapf(t,k,mk) => using(a,mk)
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk) }, { 1 + using2(a, mk) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_actions_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistency between map, mapf, labeling, and using:
    // map(t,a) & using(a,mk) => OR { mapf(t,k,mk) : k }
    void build_formulas_map_using2_then_mapf(std::ostream &os) {
        start_group("map(t,a) & using(a,mk) => OR { mapf(t,k,mk) : k }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            MetaFeature mk(tuple.at(1));
            Action a(tuple.at(2));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= a) && (a < num_actions_));

            // map(t,a) & using(a,mk) => OR { mapf(t,k,mk) : k }
            SAT::Implication IP({ 1 + map(t, a), 1 + using2(a, mk) }, { });
            for( int i = 0; i < int(features_per_layer_[layer].size()); ++i ) {
                Feature k = features_per_layer_[layer].at(i);
                IP.add_consequent(1 + mapf(t, k, mk));
            }
            add_implication(IP);
        };

        SAT::VarSet tmp(p_transitions_, p_meta_features_, p_actions_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of free(k,t,a):
    // map(t,a) & AND { -mapf(t,k,mk) : mk } => free(k,t,a)
    void build_formulas_map_mapf_then_free(std::ostream &os) {
        start_group("map(t,a) & AND { -mapf(t,k,mk) : mk } => free(k,t,a)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Transition t(tuple.at(1));
            Action a(tuple.at(2));
            Layer layer = f_layer_.at(k);
            assert((0 <= layer) && (layer < num_layers_));
            assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
            assert((0 <= a) && (a < num_actions_));

            if( layer == tr_layer_.at(t) ) {
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

                // map(t,a) & AND { -mapf(t,k,mk) : mk } => free(k,t,a)
                SAT::Implication IP({ 1 + map(t, a) }, { 1 + free(k, t, a) });
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk )
                    IP.add_antecedent(-(1 + mapf(t, k, mk)));
                add_implication(IP);
            }
        };

        free.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of free(k,t,a):
    // map(t,a) & mapf(t,k,mk) => [ -eff0(a,mk) & -eff1(a,mk) <=> free(k,t,a) ]
    void build_formulas_map_mapf_then_eff0_eff1_iff_free(std::ostream &os) {
        start_group("map(t,a) & mapf(t,k,mk) => [ -eff0(a,mk) & -eff1(a,mk) <=> free(k,t,a) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Action a(tuple.at(3));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= a) && (a < num_actions_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

            if( layer == f_layer_.at(k) ) {
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // map(t,a) & mapf(t,k,mk) & -eff0(a,mk) & -eff1(a,mk) => free(k,t,a)
                // map(t,a) & mapf(t,k,mk) & free(k,t,a) => -eff0(a,mk)
                // map(t,a) & mapf(t,k,mk) & free(k,t,a) => -eff1(a,mk)
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), -(1 + eff0(a, mk)), -(1 + eff1(a, mk)) }, { 1 + free(k, t, a) });
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + free(k, t, a) }, { -(1 + eff0(a, mk)) });
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + free(k, t, a) }, { -(1 + eff1(a, mk)) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_actions_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Transitions:
    // map(t,a) & mapf(t,k,mk) & pre0(a,mk) => -phi(k,t.src)
    // map(t,a) & mapf(t,k,mk) & pre1(a,mk) => phi(k,t.src)
    // map(t,a) & mapf(t,k,mk) & eff0(a,mk) => -phi(k,t.dst)
    // map(t,a) & mapf(t,k,mk) & eff1(a,mk) => phi(k,t.dst)
    void build_formulas_transitions(std::ostream &os) {
        start_group("(transitions) map(t,a) & mapf(t,k,mk) & pre0(a,mk) => -phi(k,t.src)");
        end_group(os, "(transitions) map(t,a) & mapf(t,k,mk) & pre0(a,mk) => -phi(k,t.src)", -1);
        start_group("(transitions) map(t,a) & mapf(t,k,mk) & eff0(a,mk) => -phi(k,t.dst)");
        end_group(os, "(transitions) map(t,a) & mapf(t,k,mk) & eff0(a,mk) => -phi(k,t.dst)", -1);
        start_group("(transitions) map(t,a) & mapf(t,k,mk) & pre1(a,mk) => phi(k,t.src)");
        end_group(os, "(transitions) map(t,a) & mapf(t,k,mk) & pre1(a,mk) => phi(k,t.src)", -1);
        start_group("(transitions) map(t,a) & mapf(t,k,mk) & eff1(a,mk) => phi(k,t.dst)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Action a(tuple.at(3));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= a) && (a < num_actions_));

            if( layer == f_layer_.at(k) ) {
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // map(t,a) & mapf(t,k,mk) & pre0(a,mk) => -phi(k,t.src)
                // map(t,a) & mapf(t,k,mk) & eff0(a,mk) => -phi(k,t.dst)
                // map(t,a) & mapf(t,k,mk) & pre1(a,mk) => phi(k,t.src)
                // map(t,a) & mapf(t,k,mk) & eff1(a,mk) => phi(k,t.dst)
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + pre0(a, mk) }, { -(1 + phi(k, tr_src_.at(t))) });
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + eff0(a, mk) }, { -(1 + phi(k, tr_dst_.at(t))) });
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + pre1(a, mk) }, { 1 + phi(k, tr_src_.at(t)) });
                add_implication({ 1 + map(t, a), 1 + mapf(t, k, mk), 1 + eff1(a, mk) }, { 1 + phi(k, tr_dst_.at(t)) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_actions_);
        tmp.enumerate_vars_from_multipliers(foo);
        int n = num_implications() - imp_offsets_.back().first;
        end_group(os, "(transitions) map(t,a) & mapf(t,k,mk) & eff1(a,mk) => phi(k,t.dst)", n / 4);
    }

    // Inertia:
    // map(t,a) => [ free(k,t,a) <=> [ phi(k,t.src) <=> phi(k,t.dst) ] ]
    void build_formulas_transitions_inertia(std::ostream &os) {
        start_group("(inertia) map(t,a) => [ free(k,t,a) <=> [ phi(k,t.src) <=> phi(k,t.dst) ] ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Transition t(tuple.at(1));
            Action a(tuple.at(2));
            Layer layer = f_layer_.at(k);
            assert((0 <= layer) && (layer < num_layers_));
            assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
            assert((0 <= a) && (a < num_actions_));

            if( layer == tr_layer_.at(t) ) {
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

                // map(t,a) & free(k,t,a) & phi(k,t.src) => phi(k,t.dst)
                // map(t,a) & free(k,t,a) & phi(k,t.dst) => phi(k,t.src)
                // map(t,a) & phi(k,t.src) & phi(k,t.dst) => free(k,t,a)
                // map(t,a) & -phi(k,t.src) & -phi(k,t.dst) => free(k,t,a)
                add_implication({ 1 + map(t, a), 1 + free(k, t, a), 1 + phi(k, tr_src_.at(t)) }, { 1 + phi(k, tr_dst_.at(t)) });
                add_implication({ 1 + map(t, a), 1 + free(k, t, a), 1 + phi(k, tr_dst_.at(t)) }, { 1 + phi(k, tr_src_.at(t)) });
                add_implication({ 1 + map(t, a), 1 + phi(k, tr_src_.at(t)), 1 + phi(k, tr_dst_.at(t)) }, { 1 + free(k, t, a) });
                add_implication({ 1 + map(t, a), -(1 + phi(k, tr_src_.at(t))), -(1 + phi(k, tr_dst_.at(t))) }, { 1 + free(k, t, a) });
            }
        };

        free.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // appl(a,t,s) => OR { mapeq(tp,a,t) : tp.src = s, tp.label = t.label }
    void build_formulas_appl_then_mapeq(std::ostream &os) {
        start_group("appl(a,t,s) => OR { mapeq(tp,a,t) : tp.src = s, tp.label = t.label }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            Transition t(tuple.at(1));
            State s(tuple.at(2));
            Layer layer = tr_layer_.at(t);
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= layer) && (layer < num_layers_));
            assert((0 <= a) && (a < num_actions_));


            if( (layer == s_layer_.at(s)) && (tr_src_.at(t) != s) ) {
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // enforce constraint only when handling a "complete" layer
                if( random_transitions_per_layer_.at(layer) == -1 ) {
                    // appl(a,t,s) => OR { mapeq(tp,a,t) : tp.src = s, tp.label = t.label }
                    SAT::Implication IP({ 1 + appl(a, t, s) }, { });
                    for( int i = 0; i < int(transitions_per_layer_[layer].size()); ++i ) {
                        Transition tp = transitions_per_layer_[layer].at(i);
                        if( (tr_src_[tp] == s) && (tr_label_[tp] == tr_label_[t]) )
                            IP.add_consequent(1 + mapeq(tp, a, t));
                    }
                    add_implication(IP);
                }
            }
        };

        appl.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // mapeq(t,a,tp) => map(t,a) & eq(t,tp)
    void build_formulas_mapeq_then_map_eq(std::ostream &os) {
        start_group("mapeq(t,a,tp) => map(t,a) & eq(t,tp)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0)), tp(tuple.at(2));
            Action a(tuple.at(1));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == tr_layer_.at(tp)) && (tr_label_.at(t) == tr_label_.at(tp)) && (t != tp) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((transitions_per_layer_.at(layer).front() <= tp) && (tp <= transitions_per_layer_.at(layer).back()));

                // mapeq(t,a,tp) => map(t,a)
                // mapeq(t,a,tp) => eq(t,tp)
                add_implication({ 1 + mapeq(t, a, tp) }, { 1 + map(t, a) });
                add_implication({ 1 + mapeq(t, a, tp) }, { 1 + eq(t < tp ? t : tp, t < tp ? tp : t) });
            }
        };

        mapeq.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // eq(t,tp) => [ mapf(t,k,mk) <=> mapf(tp,k,mk) ]
    void build_formulas_eq_then_mapf_iff_mapf(std::ostream &os) {
        start_group("eq(t,tp) => [ mapf(t,k,mk) <=> mapf(tp,k,mk) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0)), tp(tuple.at(1));
            Feature k(tuple.at(2));
            MetaFeature mk(tuple.at(3));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= mk) && (mk < num_meta_features_));

            if( (layer == tr_layer_.at(tp)) && (tr_layer_.at(t) == f_layer_.at(k)) && (tr_label_.at(t) == tr_label_.at(tp)) && (t < tp) ) {
                assert((transitions_per_layer_.at(layer).front() <= tp) && (tp <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // eq(t,tp) & mapf(t,k,mk) => mapf(tp,k,mk)
                // eq(t,tp) & mapf(tp,k,mk) => mapf(t,k,mk)
                add_implication({ 1 + eq(t, tp), 1 + mapf(t, k, mk) }, { 1 + mapf(tp, k, mk) });
                add_implication({ 1 + eq(t, tp), 1 + mapf(tp, k, mk) }, { 1 + mapf(t, k, mk) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_transitions_, p_features_, p_meta_features_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of g(k,s,t):
    // g(k,s,t) => phi(k,s) v phi(k,t)
    // g(k,s,t) => -phi(k,s) v -phi(k,t)
    // phi(k,s) & -phi(k,t) => g(k,s,t)
    // -phi(k,s) & phi(k,t) => g(k,s,t)
    void build_formulas_def_g(std::ostream &os) {
        start_group("<def g(k,s,t)>");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            State s(tuple.at(1)), t(tuple.at(2));
            if( (f_layer_.at(k) == s_layer_.at(s)) && (s_layer_.at(s) == s_layer_.at(t)) && (s < t) ) {
                Layer layer = f_layer_.at(k);
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= t) && (t <= states_per_layer_.at(layer).back()));

                // g(k,s,t) <=> [phi(k,s) & -phi(k,t)] v [-phi(k,s) & phi(k,t)]
                //          <=> [phi(k,s) v phi(k,t)] & [-phi(k,s) v -phi(k,t)]

                // g(k,s,t) => phi(k,s) v phi(k,t)
                // g(k,s,t) => -phi(k,s) v -phi(k,t)
                // phi(k,s) & -phi(k,t) => g(k,s,t)
                // -phi(k,s) & phi(k,t) => g(k,s,t)
                add_implication({ 1 + g(k, s, t) }, { 1 + phi(k, s), 1 + phi(k, t) });
                add_implication({ 1 + g(k, s, t) }, { -(1 + phi(k, s)), -(1 + phi(k, t)) });
                add_implication({ 1 + phi(k, s), -(1 + phi(k, t)) }, { 1 + g(k, s, t) });
                add_implication({ -(1 + phi(k, s)), 1 + phi(k, t) }, { 1 + g(k, s, t) });
            }
        };

        g.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Separate different states using features:
    // OR { g(k,s,t) : k } [ s < t ]
    void build_formulas_separate_states(std::ostream &os) {
        start_group("OR { g(k,s,t) : k }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            State s(tuple.at(0)), t(tuple.at(1));
            if( (s_layer_.at(s) == s_layer_.at(t)) && (s < t) ) {
                Layer layer = s_layer_.at(s);
                assert((0 <= layer) && (layer < num_layers_));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= t) && (t <= states_per_layer_.at(layer).back()));

                // OR { g(k,s,t) : k }
                SAT::Implication IP;
                for( int i = 0; i < int(features_per_layer_[layer].size()); ++i )
                    IP.add_consequent(1 + g(features_per_layer_[layer].at(i), s, t));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_states_, p_states_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of appl(a,t,s):
    // map(t,a) & AND { Z0(t,k,a,s) : k } & AND { Z1(t,k,a,s) : k } => appl(a,t,s)
    void build_formulas_map_Z0_Z1_then_appl(std::ostream &os) {
        start_group("map(t,a) & AND { Z0(t,k,a,s) & Z1(t,k,a,s) : k } => appl(a,t,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            Transition t(tuple.at(1));
            State s(tuple.at(2));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == s_layer_.at(s)) && (tr_src_.at(t) != s) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // map(t,a) & AND { Z0(t,k,a,s) : k } & AND { Z1(t,k,a,s) : k } => appl(a,t,s)
                SAT::Implication IP({ 1 + map(t, a) }, { 1 + appl(a, t, s) });
                for( int i = 0; i < int(features_per_layer_[layer].size()); ++i ) {
                    IP.add_antecedent(1 + Z0(t, features_per_layer_[layer].at(i), a, s));
                    IP.add_antecedent(1 + Z1(t, features_per_layer_[layer].at(i), a, s));
                }
                add_implication(IP);
            }
        };

        appl.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z0(t,k,a,s): [ mapf0(t,k,a) => -phi(k,s) ] => Z0(t,k,a,s)
    // -phi(k,s) => Z0(t,k,a,s)
    void build_formulas_phi_then_Z0(std::ostream &os) {
        start_group("-phi(k,s) => Z0(t,k,a,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Action a(tuple.at(2));
            State s(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (f_layer_.at(k) == s_layer_.at(s)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // -phi(k,s) => Z0(t,k,a,s)
                add_implication({ -(1 + phi(k, s)) }, { 1 + Z0(t, k, a, s) });
            }
        };

        Z0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z0(t,k,a,s): [ mapf0(t,k,a) => -phi(k,s) ] => Z0(t,k,a,s)
    // AND { -X0(a,t,k,mk) : mk } => Z0(t,k,a,s)
    void build_formulas_X0_then_Z0(std::ostream &os) {
        start_group("AND { -X0(a,t,k,mk) : mk } => Z0(t,k,a,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Action a(tuple.at(2));
            State s(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (f_layer_.at(k) == s_layer_.at(s)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // AND { -X0(a,t,k,mk) : mk } => Z0(t,k,a,s)
                SAT::Implication IP({ }, { 1 + Z0(t, k, a, s) });
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk )
                    IP.add_antecedent(-(1 + X0(a, t, k, mk)));
                add_implication(IP);
            }
        };

        Z0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z0(t,k,a,s): [ mapf0(t,k,a) => -phi(k,s) ] => Z0(t,k,a,s)
    // X0(a,t,k,mk) => pre0(a,mk) & mapf(t,k,mk)
    void build_formulas_X0_then_pre0_mapf(std::ostream &os) {
        start_group("X0(a,t,k,mk) => pre0(a,mk) & mapf(t,k,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            Transition t(tuple.at(1));
            Feature k(tuple.at(2));
            MetaFeature mk(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            if( tr_layer_.at(t) == f_layer_.at(k) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // X0(a,t,k,mk) => pre0(a,mk)
                // X0(a,t,k,mk) => mapf(t,k,mk)
                add_implication({ 1 + X0(a, t, k, mk) }, { 1 + pre0(a, mk) });
                add_implication({ 1 + X0(a, t, k, mk) }, { 1 + mapf(t, k, mk) });
            }
        };

        X0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z1(t,k,a,s): [ mapf1(t,k,a) => phi(k,s) ] => Z1(t,k,a,s)
    // phi(k,s) => Z1(t,k,a,s)
    void build_formulas_phi_then_Z1(std::ostream &os) {
        start_group("phi(k,s) => Z1(t,k,a,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Action a(tuple.at(2));
            State s(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (f_layer_.at(k) == s_layer_.at(s)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // phi(k,s) => Z1(t,k,a,s)
                add_implication({ 1 + phi(k, s) }, { 1 + Z1(t, k, a, s) });
            }
        };

        Z1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z1(t,k,a,s): [ mapf1(t,k,a) => phi(k,s) ] => Z1(t,k,a,s)
    // AND { -X1(a,t,k,mk) : mk } => Z1(t,k,a,s)
    void build_formulas_X1_then_Z1(std::ostream &os) {
        start_group("AND { -X1(a,t,k,mk) : mk } => Z1(t,k,a,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Action a(tuple.at(2));
            State s(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (f_layer_.at(k) == s_layer_.at(s)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // AND { -X1(a,t,k,mk) : mk } => Z1(t,k,a,s)
                SAT::Implication IP({ }, { 1 + Z1(t, k, a, s) });
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk )
                    IP.add_antecedent(-(1 + X1(a, t, k, mk)));
                add_implication(IP);
            }
        };

        Z1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of Z1(t,k,a,s): [ mapf1(t,k,a) => phi(k,s) ] => Z1(t,k,a,s)
    // X1(a,t,k,mk) => pre1(a,mk) & mapf(t,k,mk)
    void build_formulas_X1_then_pre1_mapf(std::ostream &os) {
        start_group("X1(a,t,k,mk) => pre1(a,mk) & mapf(t,k,mk)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            Transition t(tuple.at(1));
            Feature k(tuple.at(2));
            MetaFeature mk(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mk) && (mk < num_meta_features_));

            if( tr_layer_.at(t) == f_layer_.at(k) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // X1(a,t,k,mk) => pre1(a,mk)
                // X1(a,t,k,mk) => mapf(t,k,mk)
                add_implication({ 1 + X1(a, t, k, mk) }, { 1 + pre1(a, mk) });
                add_implication({ 1 + X1(a, t, k, mk) }, { 1 + mapf(t, k, mk) });
            }
        };

        X1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // EXACTLY-1 { ground(k,p) : p }
    void build_formulas_exactly_1_ground1(std::ostream &os) {
        start_group("EXACTLY-1 { ground(k,p) : p }");

        for( Feature k = Feature(0); k < num_features_; ++k ) {
            std::vector<int> literals(num_atoms_);
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                literals[p] = 1 + ground1(k, p);
            exactly_1(std::string("ground-exactly-1(") + std::to_string(k) + ")", literals);
        }

        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // AT-MOST-1 { ground(k,i,o) : o }
    void build_formulas_at_most_1_ground2(std::ostream &os) {
        start_group("AT-MOST-1 { ground(k,i,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Arity i(tuple.at(1));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( 0 < i ) {
                Layer layer = f_layer_.at(k);
                assert((0 <= layer) && (layer < num_layers_));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                std::vector<int> literals(objects_per_layer_.at(layer));
                for( Object o = Object(0); o < objects_per_layer_.at(layer); ++o )
                    literals[o] = 1 + ground2(k, i, o);
                at_most_1(std::string("ground-at-most-1(") + std::to_string(k) + "," + std::to_string(i) + ")", literals);
            }
        };

        SAT::VarSet tmp(p_features_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // ground(k,p) & ground(k,i,o) => OR { arity(p,j) : i <= j <= max-arity }
    void build_formulas_ground1_ground2_then_arity(std::ostream &os) {
        start_group("ground(k,p) & ground(k,i,o) => OR { arity(p,j) : i <= j }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2));
            Object o(tuple.at(3));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (0 < i) && (o < objects_per_layer_.at(f_layer_.at(k))) ) {
                Layer layer = f_layer_.at(k);
                assert((0 <= layer) && (layer < num_layers_));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));
                assert((0 <= o) && (o < objects_per_layer_.at(layer)));

                // ground(mk,p) & ground(mk,i,mo) => OR { arity(p,j) : i <= j <= max-arity }
                SAT::Implication IP({ 1 + ground1(k, p), 1 + ground2(k, i, o) }, { });
                for( Arity j = i; j < 1 + max_arity_; ++j )
                    IP.add_consequent(1 + arity(p, j));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_features_, p_atoms_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // ground(k,p) & arity(p,i) => OR { ground(k,j,o) : o } [ 1 <= j <= i ]
    void build_formulas_ground1_arity_then_ground2(std::ostream &os) {
        start_group("ground(k,p) & arity(p,i) => OR { ground(k,j,o) : o } [1 <= j <= i ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2)), j(tuple.at(3));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (0 < j) && (j <= i) ) {
                Layer layer = f_layer_.at(k);
                assert((0 <= layer) && (layer < num_layers_));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // ground(k,p) & arity(p,i) => OR { ground(k,j,o) : o } [ 1 <= j <= i ]
                SAT::Implication IP({ 1 + ground1(k, p), 1 + arity(p, i) }, { });
                for( Object o = Object(0); o < objects_per_layer_.at(layer); ++o )
                    IP.add_consequent(1 + ground2(k, j, o));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_features_, p_atoms_, p_arities_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // ground(k,p) & arity(p,i) => -ground(k,j,o) [ i < j ]
    void build_formulas_ground1_arity_then_ground2_2(std::ostream &os) {
        start_group("ground(k,p) & arity(p,i) => -ground(k,j,o) [ i < j ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Atom p(tuple.at(1));
            Arity i(tuple.at(2)), j(tuple.at(3));
            Object o(tuple.at(4));
            assert((0 <= k) && (k < num_features_));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= j) && (j < 1 + max_arity_));

            if( (i < j) && (o < objects_per_layer_.at(f_layer_.at(k))) ) {
                // ground(k,p) & arity(p,i) => -ground(k,j,o) [ i < j ]
                add_implication({ 1 + ground1(k, p), 1 + arity(p, i) }, { -(1 + ground2(k, j, o)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_atoms_, p_arities_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Mapping of features to grounded atoms:
    // -ground(k,0,o)
    void build_formulas_ground2(std::ostream &os) {
        start_group("-ground(k,0,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0));
            Object o(tuple.at(1));
            assert((0 <= k) && (k < num_features_));

            if( o < objects_per_layer_.at(f_layer_.at(k)) ) {
                // -ground(k,0,o)
                add_unit(-(1 + ground2(k, 0, o)));
            }
        };

        SAT::VarSet tmp(p_features_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc1) ground(k,p) & ground(kp,p) & arity(p,i) => OR { gdiff(k,kp,j) : 1 <= j <= i }
    void ENC1_build_formulas_ground1_ground1_arity_then_gdiff(std::ostream &os) {
        start_group("(enc1) ground(k,p) & ground(kp,p) & arity(p,i) => OR { gdiff(k,kp,j) : 1 <= j <= i }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            Arity i(tuple.at(3));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // ground(k,p) & ground(kp,p) & arity(p,i) => OR { gdiff(k,kp,j) : 1 <= j <= i }
                SAT::Implication IP({ 1 + ground1(k, p), 1 + ground1(kp, p), 1 + arity(p, i) }, { });
                for( Arity j = Arity(1); j <= i; ++j )
                    IP.add_consequent(1 + gdiff(k, kp, j));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_atoms_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc1) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
    void ENC1_build_formulas_gdiff_then_ground2_ground2(std::ostream &os) {
        start_group("(enc1) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Arity i(tuple.at(2));
            Object o(tuple.at(3));
            assert((0 <= i) && (i < 1 + max_arity_));
            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) && (o < objects_per_layer_.at(f_layer_.at(k))) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((0 <= o) && (o < objects_per_layer_.at(f_layer_.at(k))));

                // gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
                add_implication({ 1 + gdiff(k, kp, i) }, { -(1 + ground2(k, i, o)), -(1 + ground2(kp, i, o)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc1) gdiff(k,kp,i) => OR { ground(k,i,o) : o } & OR { ground(kp,i,o) : o }
    void ENC1_build_formulas_gdiff_then_ground2_ground2_2(std::ostream &os) {
        start_group("(enc1) gdiff(k,kp,i) => OR { ground(k,i,o) : o } & OR { ground(kp,i,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Arity i(tuple.at(2));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,i) => OR { ground(k,i,o) : o }
                SAT::Implication IP1({ 1 + gdiff(k, kp, i) }, { });
                for( Object o = Object(0); o < objects_per_layer_.at(f_layer_.at(k)); ++o )
                    IP1.add_consequent(1 + ground2(k, i, o));
                add_implication(IP1);

                // gdiff(k,kp,i) => OR { ground(kp,i,o) : o }
                SAT::Implication IP2({ 1 + gdiff(k, kp, i) }, { });
                for( Object o = Object(0); o < objects_per_layer_.at(f_layer_.at(k)); ++o )
                    IP2.add_consequent(1 + ground2(kp, i, o));
                add_implication(IP2);
            }
        };

        gdiff.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc1) gdiff(k,kp,i) => AND { -gdiff(k,kp,j) : 1 <= j < i }
    void ENC1_build_formulas_gdiff_then_gdiff(std::ostream &os) {
        start_group("(enc1) gdiff(k,kp,i) => AND { -gdiff(k,kp,j) : 1 <= j < i }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Arity i(tuple.at(2));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,i) => AND { -gdiff(k,kp,j) : 1 <= j < i }
                for( Arity j = Arity(1); j < i; ++j )
                    add_implication({ 1 + gdiff(k, kp, i) }, { -(1 + gdiff(k, kp, j)) });
            }
        };

        gdiff.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc2) arity(p,0) => -ground(k,p) v -ground(kp,p)
    void ENC2_build_formulas_arity_then_ground1_ground1(std::ostream &os) {
        start_group("(enc2) arity(p,0) => -ground(k,p) v -ground(kp,p)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            assert((0 <= p) && (p < num_atoms_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // arity(p,0) => -ground(k,p) v -ground(kp,p)
                add_implication({ 1 + arity(p, 0) }, { -(1 + ground1(k, p)), -(1 + ground1(kp, p)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_atoms_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc2) OR { gdiff(k,kp,i) : 0 <= i <= max-arity }
    void ENC2_build_formulas_gdiff(std::ostream &os) {
        start_group("(enc2) OR { gdiff(k,kp,i) : 0 <= i <= max-arity }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // OR { gdiff(k,kp,i) : 0 <= i <= max-arity }
                SAT::Implication IP;
                for( Arity i = Arity(0); i < 1 + max_arity_; ++i )
                    IP.add_consequent(1 + gdiff(k, kp, i));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_features_, p_features_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc2) gdiff(k,kp,0) => -ground(k,p) v -ground(kp,p)
    void ENC2_build_formulas_gdiff_then_ground1_ground1(std::ostream &os) {
        start_group("(enc2) gdiff(k,kp,0) => -ground(k,p) v -ground(kp,p)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            assert((0 <= p) && (p < num_atoms_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,0) => -ground(k,p) v -ground(kp,p)
                add_implication({ 1 + gdiff(k, kp, 0) }, { -(1 + ground1(k, p)), -(1 + ground1(kp, p)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_atoms_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc2) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
    void ENC2_build_formulas_gdiff_then_ground2_ground2(std::ostream &os) {
        start_group("(enc2) gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Arity i(tuple.at(2));
            Object o(tuple.at(3));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) && (o < objects_per_layer_.at(f_layer_.at(k))) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((0 <= o) && (o < objects_per_layer_.at(f_layer_.at(k))));

                // gdiff(k,kp,i) => -ground(k,i,o) v -ground(kp,i,o)
                add_implication({ 1 + gdiff(k, kp, i) }, { -(1 + ground2(k, i, o)), -(1 + ground2(kp, i, o)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc2) gdiff(k,kp,i) => OR { ground(k,i,o) : o } v OR { ground(kp,i,o) : o }
    void ENC2_build_formulas_gdiff_then_ground2_ground2_2(std::ostream &os) {
        start_group("(enc2) gdiff(k,kp,i) => OR { ground(k,i,o) : o } v OR { ground(kp,i,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Arity i(tuple.at(2));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,i) => OR { ground(k,i,o) : o } v OR { ground(kp,i,o) : o }
                SAT::Implication IP({ 1 + gdiff(k, kp, i) }, { });
                for( Object o = Object(0); o < objects_per_layer_.at(f_layer_.at(k)); ++o ) {
                    IP.add_consequent(1 + ground2(k, i, o));
                    IP.add_consequent(1 + ground2(kp, i, o));
                }
                add_implication(IP);
            }
        };

        gdiff.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (enc2) [DISABLED] AT-MOST-1 { gdiff(k,kp,i) : 0 <= i <= max-arity }
    void ENC2_build_formulas_at_most_1_gdiff(std::ostream &os) {
        start_group("(enc2) AT-MOST-1 { gdiff(k,kp,i) : 0 <= i <= max-arity }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // AT-MOST-1 { gdiff(k,kp,i) : 0 <= i <= max-arity }
                std::vector<int> literals(1 + max_arity_);
                for( Arity i = Arity(0); i < 1 + max_arity_; ++i )
                    literals[i] = 1 + gdiff(k, kp, i);
                at_most_1(std::string("at-most-1-gdiff(") + std::to_string(k) + "," + std::to_string(kp) + ")", literals);
            }
        };

        SAT::VarSet tmp(p_features_, p_features_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc3) gdiff(k,kp,p) <=> ground(k,p) & ground(kp,p)
    void ENC3_build_formulas_gdiff1_iff_ground1_ground1(std::ostream &os) {
        start_group("(enc3) gdiff(k,kp,p) <=> ground(k,p) & ground(kp,p)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            assert((0 <= p) && (p < num_atoms_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,p) => ground(k,p)
                // gdiff(k,kp,p) => ground(kp,p)
                // ground(k,p) & ground(kp,p) => gdiff(k,kp,p)
                add_implication({ 1 + gdiff1(k, kp, p) }, { 1 + ground1(k, p) });
                add_implication({ 1 + gdiff1(k, kp, p) }, { 1 + ground1(kp, p) });
                add_implication({ 1 + ground1(k, p), 1 + ground1(kp, p) }, { 1 + gdiff1(k, kp, p) });
            }
        };

        gdiff1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc3) gdiff(k,kp,p) & arity(p,i) => OR { gdiff(k,kp,p,j) : 1 <= j <= i }
    void ENC3_build_formulas_gdiff1_arity_then_gdiff2(std::ostream &os) {
        start_group("(enc3) gdiff(k,kp,p) & arity(p,i) => OR { gdiff(k,kp,p,j) : 1 <= j <= i }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            Arity i(tuple.at(3));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,p) & arity(p,i) => OR { gdiff(k,kp,p,j) : 1 <= j <= i }
                SAT::Implication IP({ 1 + gdiff1(k, kp, p), 1 + arity(p, i) }, { });
                for( Arity j = Arity(1); j <= i; ++j )
                    IP.add_consequent(1 + gdiff2(k, kp, p, j));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_atoms_, p_arities_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc3) gdiff(k,kp,p,j) => gdiff(k,kp,p) & OR { arity(p,i) : j <= i <= max-arity }
    void ENC3_build_formulas_gdiff2_then_gdiff1_arity(std::ostream &os) {
        start_group("(enc3) gdiff(k,kp,p,i,j) => gdiff(k,kp,p) & OR { arity(p,i) : j <= i <= max-arity }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            Arity j(tuple.at(3));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= j) && (j < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < j) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));

                // gdiff(k,kp,p,j) => gdiff(k,kp,p)
                add_implication({ 1 + gdiff2(k, kp, p, j) }, { 1 + gdiff1(k, kp, p) });

                // gdiff(k,kp,p,j) => OR { arity(p,i) : j <= i <= max-arity }
                SAT::Implication IP({ 1 + gdiff2(k, kp, p, j) }, { });
                for( Arity i = Arity(j); i < 1 + max_arity_; ++i )
                    IP.add_consequent(1 + arity(p, i));
                add_implication(IP);
            }
        };

        gdiff2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to ground atoms:
    // (enc3) gdiff(k,kp,p,i) => -ground(k,i,o) v -ground(kp,i,o)
    void ENC3_build_formulas_gdiff2_then_ground2_ground2(std::ostream &os) {
        start_group("(enc3) gdiff(k,kp,p,i) => -ground(k,i,o) v -ground(kp,i,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Feature k(tuple.at(0)), kp(tuple.at(1));
            Atom p(tuple.at(2));
            Arity i(tuple.at(3));
            Object o(tuple.at(4));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (f_layer_.at(k) == f_layer_.at(kp)) && (k < kp) && (0 < i) && (o < objects_per_layer_.at(f_layer_.at(k))) ) {
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= k) && (k <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((features_per_layer_.at(f_layer_.at(k)).front() <= kp) && (kp <= features_per_layer_.at(f_layer_.at(k)).back()));
                assert((0 <= o) && (o < objects_per_layer_.at(f_layer_.at(k))));

                // gdiff(k,kp,p,i) => -ground(k,i,o) v -ground(kp,i,o)
                add_implication({ 1 + gdiff2(k, kp, p, i) }, { -(1 + ground2(k, i, o)), -(1 + ground2(kp, i, o)) });
            }
        };

        SAT::VarSet tmp(p_features_, p_features_, p_atoms_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // 1-1 mapping of features to grounded atoms:
    // (enc123) <non-strict-lexicographic-ordering-features-in-layers>
    // (enc4) <strict-lexicographic-ordering-features-in-layers>
    void build_formulas_lex_ordering_features_in_layers(std::ostream &os) {
        assert(some_features_to_actions_map_encoding());
        if( !some_features_to_actions_map_encoding() )
            start_group("<lexicographic-ordering-features-in-layers>");
        else if( !encoding("features-to-actions-map-strict-lex-ordering") )
            start_group("(enc123) <non-strict-lexicographic-ordering-features-in-layers>");
        else if( encoding("features-to-actions-map-strict-lex-ordering") )
            start_group("(enc4) <strict-lexicographic-ordering-features-in-layers>");

        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            std::string prefix(std::string("features-in-layer-") + std::to_string(layer));
            std::vector<std::vector<int> > lvectors;

            for( int i = 0; i < int(features_per_layer_.at(layer).size()); ++i ) {
                Feature k = Feature(features_per_layer_.at(layer).at(i));
                assert(f_layer_.at(k) == layer);
                std::vector<int> feature;

                if( encoding("features-to-actions-map-strict-lex-ordering") ) {
                    // order by grounded atom names
                    for( Atom p = Atom(0); p < num_atoms_; ++p )
                        feature.push_back(-(1 + ground1(k, p)));

                    // order by arguments
                    for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                        for( Object o = Object(0); o < objects_per_layer_.at(layer); ++o )
                            feature.push_back(-(1 + ground2(k, i, o)));
                    }
                } else {
                    for( int j = 0; j < int(states_per_layer_.at(layer).size()); ++j ) {
                        State s = State(states_per_layer_.at(layer).at(j));
                        feature.push_back(1 + phi(k, s));
                    }
                }

                lvectors.emplace_back(std::move(feature));
            }

            if( !lvectors.empty() ) lex_ordering(prefix, lvectors, encoding("features-to-actions-map-strict-lex-ordering"));
        }

        end_group(os);
    }

    // Consistent mapping between features and meta-features
    // mapf(t,k,mk) => [ atom(mk,p) <=> ground(k,p) ]
    void build_formulas_mapf_then_atom1_iff_ground1(std::ostream &os) {
        start_group("mapf(t,k,mk) => [ atom(mk,p) <=> ground(k,p) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            assert((0 <= mk) && (mk < num_meta_features_));

            if( tr_layer_.at(t) == f_layer_.at(k) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // mapf(t,k,mk) => [ atom(mk,p) <=> ground(k,p) ]
                for( Atom p = Atom(0); p < num_atoms_; ++p ) {
                    // mapf(t,k,mk) & atom(mk,p) => ground(k,p)
                    // mapf(t,k,mk) & ground(k,p) => atom(mk,p)
                    add_implication({ 1 + mapf(t, k, mk), 1 + atom1(mk, p) }, { 1 + ground1(k, p) });
                    add_implication({ 1 + mapf(t, k, mk), 1 + ground1(k, p) }, { 1 + atom1(mk, p) });
                }
            }
        };

        mapf.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent mapping between features and meta-features
    // mapf(t,k,mk) & atom(mk,i,mo) => OR { ground(k,i,o) : o }
    void build_formulas_mapf_atom2_then_ground2(std::ostream &os) {
        start_group("mapf(t,k,mk) & atom(mk,i,mo) => OR { ground(k,i,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            MetaObject mo(tuple.at(4));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (0 < i) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // mapf(t,k,mk) & atom(mk,i,mo) => OR { ground(k,i,o) : o }
                SAT::Implication IP({ 1 + mapf(t, k, mk), 1 + atom2(mk, i, mo) }, { });
                for( Object o = Object(0); o < objects_per_layer_.at(tr_layer_.at(t)); ++o )
                    IP.add_consequent(1 + ground2(k, i, o));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_arities_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Consistent mapping between features and meta-features
    // mapf(t,k,mk) & ground(k,i,o) => OR { atom(km,i,mo) : mo }
    void build_formulas_mapf_ground2_then_atom2(std::ostream &os) {
        start_group("mapf(t,k,mk) & ground(k,i,o) => OR { atom(mk,i,mo) : mo }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            Object o(tuple.at(4));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (0 < i) && (o < objects_per_layer_.at(tr_layer_.at(t))) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // mapf(t,k,mk) & ground(k,i,o) => OR { atom(km,i,mo) : mo }
                SAT::Implication IP({ 1 + mapf(t, k, mk), 1 + ground2(k, i, o) }, { });
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    IP.add_consequent(1 + atom2(mk, i, mo));
                add_implication(IP);
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_arities_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of U(l,u,a,mo,o) and B(l,b,a,mo,mop,o,op):
    // U(l,u,a,mo,o) <=> unary(u,a,mo) & -r(l,u,o)
    void build_formulas_U_iff_unary_r(std::ostream &os) {
        start_group("U(l,u,a,mo,o) <=> unary(u,a,mo) & -r(l,u,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Layer l(tuple.at(0));
            Unary u(tuple.at(1));
            Action a(tuple.at(2));
            MetaObject mo(tuple.at(3));
            Object o(tuple.at(4));
            assert((0 <= l) && (l < num_layers_));
            assert((0 <= u) && (u < num_static_unary_predicates_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( o < objects_per_layer_.at(l) ) {
                // U(l,u,a,mo,o) => unary(u,a,mo)
                // U(l,u,a,mo,o) => -r(l,u,o)
                // U(l,u,a,mo,o) <= unary(u,a,mo) & -r(l,u,o)
                add_implication({ 1 + U(l, u, a, mo, o) }, { 1 + unary(u, a, mo) });
                add_implication({ 1 + U(l, u, a, mo, o) }, { -(1 + r(l, u, o)) });
                add_implication({ 1 + unary(u, a, mo), -(1 + r(l, u, o)) }, { 1 + U(l, u, a, mo, o) });
            }
        };

        U.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of U(l,u,a,mo,o) and B(l,b,a,mo,mop,o,op):
    // B(l,b,a,mo,mop,o,op) <=> binary(b,a,mo,mop) & -s(l,b,o,op)
    void build_formulas_B_iff_binary_s(std::ostream &os) {
        start_group("B(l,b,a,mo,mop,o,op) <=> binary(b,a,mo,mop) & -s(l,b,o,op)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Layer l(tuple.at(0));
            Binary b(tuple.at(1));
            Action a(tuple.at(2));
            MetaObject mo(tuple.at(3)), mop(tuple.at(4));
            Object o(tuple.at(5)), op(tuple.at(6));
            assert((0 <= l) && (l < num_layers_));
            assert((0 <= b) && (b < num_static_binary_predicates_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= mop) && (mop < num_meta_objects_));

            if( (mo < mop) && (o < objects_per_layer_.at(l)) && (op < objects_per_layer_.at(l)) ) {
                // B(l,b,a,mo,mop,o,op) => binary(b,a,mo,mop)
                // B(l,b,a,mo,mop,o,op) => -s(l,b,o,op)
                // B(l,b,a,mo,mop,o,op) <= binary(b,a,mo,mop) & -s(l,b,o,op)
                add_implication({ 1 + B(l, b, a, mo, mop, o, op) }, { 1 + binary(b, a, mo, mop) });
                add_implication({ 1 + B(l, b, a, mo, mop, o, op) }, { -(1 + s(l, b, o, op)) });
                add_implication({ 1 + binary(b, a, mo, mop), -(1 + s(l, b, o, op)) }, { 1 + B(l, b, a, mo, mop, o, op) });
            }
        };

        B.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of mapt(t,mo,o):
    // AT-MOST-1 { mapt(t,mo,o) : o }
    void build_formulas_at_most_1_mapt(std::ostream &os) {
        start_group("AT-MOST-1 { mapt(t,mo,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            MetaObject mo(tuple.at(1));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
            assert((0 <= mo) && (mo < num_meta_objects_));

            // AT-MOST-1 { mapt(t,mo,o) : o }
            std::vector<int> literals(objects_per_layer_[layer]);
            for( Object o = Object(0); o < objects_per_layer_[layer]; ++o )
                literals[o] = 1 + mapt(t, mo, o);
            at_most_1(std::string("mapt-at-most-1(") + std::to_string(t) + "," + std::to_string(mo) + ")", literals);
        };

        SAT::VarSet tmp(p_transitions_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of mapt(t,mo,o):
    // map(t,a) & args(a,mo) => OR { mapt(t,mo,o) : o }
    void build_formulas_map_args_then_mapt(std::ostream &os) {
        start_group("map(t,a) & args(a,mo) => OR { mapt(t,mo,o) : o }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2));
            Layer layer = tr_layer_.at(t);
            assert((0 <= layer) && (layer < num_layers_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

            // map(t,a) & args(a,mo) => OR { mapt(t,mo,o) : o }
            SAT::Implication IP({ 1 + map(t, a), 1 + args(a, mo) }, { });
            for( Object o = Object(0); o < objects_per_layer_[layer]; ++o )
                IP.add_consequent(1 + mapt(t, mo, o));
            add_implication(IP);
        };

        SAT::VarSet tmp(p_transitions_, p_actions_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of mapt(t,mo,o):
    // map(t,a) & mapt(t,mo,o) => args(a,mo)
    void build_formulas_map_mapt_then_args(std::ostream &os) {
        start_group("map(t,a) & mapt(t,mo,o) => args(a,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2));
            Object o(tuple.at(3));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( o < objects_per_layer_.at(tr_layer_.at(t)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((0 <= 0) && (o < objects_per_layer_.at(layer)));

                //  map(t,a) & mapt(t,mo,o) => args(a,mo)
                add_implication({ 1 + map(t, a), 1 + mapt(t, mo, o) }, { 1 + args(a, mo) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_actions_, p_meta_objects_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Cross consistency between schemas and transitions:
    // mapf(t,k,mk) & atom(mk,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]
    void build_formulas_mapf_atom2_then_ground2_iff_mapt(std::ostream &os) {
        start_group("mapf(t,k,mk) & atom(mk,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            MetaObject mo(tuple.at(4));
            Object o(tuple.at(5));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (0 < i) && (o < objects_per_layer_.at(tr_layer_.at(t))) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // mapf(t,k,mk) & atom(mk,i,mo) & ground(k,i,o) => mapt(t,mo,o)
                // mapf(t,k,mk) & atom(mk,i,mo) & mapt(t,mo,o) => ground(k,i,o)
                add_implication({ 1 + mapf(t, k, mk), 1 + atom2(mk, i, mo), 1 + ground2(k, i, o) }, { 1 + mapt(t, mo, o) });
                add_implication({ 1 + mapf(t, k, mk), 1 + atom2(mk, i, mo), 1 + mapt(t, mo, o) }, { 1 + ground2(k, i, o) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_arities_, p_meta_objects_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Cross consistency between schemas and transitions:
    // mapf(t,k,mk) & atom(mk,i,mo) => W(t,k,i,mo)
    void build_formulas_mapf_atom2_then_W(std::ostream &os) {
        start_group("mapf(t,k,mk) & atom(mk,i,mo) => W(t,k,i,mo)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            MetaFeature mk(tuple.at(2));
            Arity i(tuple.at(3));
            MetaObject mo(tuple.at(4));
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (0 < i) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // mapf(t,k,mk) & atom(mk,i,mo) => W(t,k,i,mo)
                add_implication({ 1 + mapf(t, k, mk), 1 + atom2(mk, i, mo) }, { 1 + W(t, k, i, mo) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_meta_features_, p_arities_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Cross consistency between schemas and transitions:
    // W(t,k,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]
    void build_formulas_W_then_ground2_iff_mapt(std::ostream &os) {
        start_group("W(t,k,i,mo) => [ ground(k,i,o) <=> mapt(t,mo,o) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Feature k(tuple.at(1));
            Arity i(tuple.at(2));
            MetaObject mo(tuple.at(3));
            Object o(tuple.at(4));
            assert((0 <= i) && (i < 1 + max_arity_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( (tr_layer_.at(t) == f_layer_.at(k)) && (0 < i) && (o < objects_per_layer_.at(tr_layer_.at(t))) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));
                assert((features_per_layer_.at(layer).front() <= k) && (k <= features_per_layer_.at(layer).back()));

                // W(t,k,i,mo) & ground(k,i,o) => mapt(t,mo,o)
                // W(t,k,i,mo) & mapt(t,mo,o) => ground(k,i,o)
                add_implication({ 1 + W(t, k, i, mo), 1 + ground2(k, i, o) }, { 1 + mapt(t, mo, o) });
                add_implication({ 1 + W(t, k, i, mo), 1 + mapt(t, mo, o) }, { 1 + ground2(k, i, o) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_features_, p_arities_, p_meta_objects_, p_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // One-hot encoding of variables:
    // AT-LEAST-<LB> { phi(k,s) : k } [ for s in layer where <LB> is LB on feature sum per state (if != 0) ]
    // AT-MOST-<UB> { phi(k,s) : k } [ for s in layer where <UB> is UB on feature sum per state (if != 0) ]
    void build_formulas_feature_sum(std::ostream &os) {
        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            start_group(std::string("<bounds-feature-sum-layer-") + std::to_string(layer) + ">");

            for( int i = 0; i < int(states_per_layer_[layer].size()); ++i ) {
                State s = states_per_layer_[layer].at(i);
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                int num_features_in_layer = features_per_layer_.at(layer).size();
                std::vector<int> literals(num_features_in_layer);
                for( int i = 0; i < num_features_in_layer; ++i ) {
                    Feature k = features_per_layer_.at(layer).at(i);
                    literals[i] = 1 + phi(k, s);
                }

                int lower_bound = feature_sum_per_layer_lower_[layer];
                int upper_bound = feature_sum_per_layer_upper_[layer];
                if( (lower_bound > 0) || (upper_bound > 0) ) {
                    std::vector<int> z;
                    sorting_network(std::string("feature-sum(") + std::to_string(s) + ")", literals, z);
                    assert(z.size() == literals.size());

                    // enforce lower bound (if any)
                    if( lower_bound > 0 ) {
                        if( lower_bound > int(literals.size()) ) {
                            add_empty_clause();
                        } else {
                            for( int i = 0; i < lower_bound; ++i )
                                add_unit(1 + z.at(i));
                        }
                    }

                    // enforce upper bound (if any)
                    if( upper_bound > 0 ) {
                        if( upper_bound < int(literals.size()) ) {
                            for( int i = upper_bound; i < int(literals.size()); ++i )
                                add_unit(-(1 + z.at(i)));
                        }
                    }
                }
            }

            end_group(os);
        }
    }

    // One-hot encoding of variables:
    // (enc123) -phi(k,s0) [ k < #features - <upper-bound> ]
    void build_formulas_negative_units_features(std::ostream &os) {
        assert(some_features_to_actions_map_encoding() && !encoding("features-to-actions-map-strict-lex-ordering"));
        start_group("-phi(k,s0) [ k < #features - <upper-bound> ]");

        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            int num_features_in_layer = features_per_layer_.at(layer).size();
            int upper_bound = feature_sum_per_layer_upper_.at(layer);
            assert(upper_bound >= 0);
            if( upper_bound > 0 ) {
                for( int i = 0; i < num_features_in_layer - upper_bound; ++i ) {
                    if( i < int(features_per_layer_.at(layer).size()) ) {
                        Feature k = features_per_layer_.at(layer).at(i);
                        assert((0 <= k) && (k < num_features_));
                        assert(k < num_features_ - feature_sum_per_layer_upper_.at(f_layer_.at(k)));

                        int s0 = states_per_layer_.at(f_layer_.at(k)).front();
                        add_unit(-(1 + phi(k, s0)));
                    }
                }
            }
        }

        end_group(os);
    }

    // One-hot encoding of variables:
    // (enc123) phi(k,s0) [ k >= #features - <lower-bound> ]
    void build_formulas_positive_units_features(std::ostream &os) {
        assert(some_features_to_actions_map_encoding() && !encoding("features-to-actions-map-strict-lex-ordering"));
        start_group("phi(k,s0) [ k >= #features - <lower-bound> ]");

        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            int num_features_in_layer = features_per_layer_.at(layer).size();
            int lower_bound = feature_sum_per_layer_lower_.at(layer);
            assert(lower_bound >= 0);
            if( lower_bound > 0 ) {
                for( int i = std::max(num_features_in_layer - lower_bound, 0); i < num_features_in_layer; ++i ) {
                    Feature k = features_per_layer_.at(layer).at(i);
                    assert((0 <= k) && (k < num_features_));
                    assert(k >= std::max(0, int(features_per_layer_.at(f_layer_.at(k)).size()) - feature_sum_per_layer_lower_.at(f_layer_.at(k))));

                    int s0 = states_per_layer_.at(f_layer_.at(k)).front();
                    add_unit(1 + phi(k, s0));

                }
            }
        }

        end_group(os);
    }

    // Explanation of existing grounded actions (common to different encodings):
    // map(t,a) & mapt(t,mo,o) & unary(u,a,mo) => r(t.layer,u,o)
    void build_formulas_map_mapt_unary_then_r(std::ostream &os) {
        start_group("map(t,a) & mapt(t,mo,o) & unary(u,a,mo) => r(t.layer,u,o)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2));
            Object o(tuple.at(3));
            Unary u(tuple.at(4));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= u) && (u < num_static_unary_predicates_));

            if( o < objects_per_layer_.at(tr_layer_.at(t)) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

                // map(t,a) & mapt(t,mo,o) & unary(u,a,mo) => r(t.layer,u,o)
                add_implication({ 1 + map(t, a), 1 + mapt(t, mo, o), 1 + unary(u, a, mo) }, { 1 + r(layer, u, o) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_actions_, p_meta_objects_, p_objects_, p_unary_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Explanation of existing grounded actions (common to different encodings):
    // map(t,a) & mapt(t,mo,o) & mapt(t,mop,op) & binary(b,a,mo,mop) => s(t.layer,b,o,op)
    void build_formulas_map_mapt_mapt_binary_then_s(std::ostream &os) {
        start_group("map(t,a) & mapt(t,mo,o) & mapt(t,mop,op) & binary(t.layer,b,a,mo,mop) => s(t.layer,b,o,op)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Transition t(tuple.at(0));
            Action a(tuple.at(1));
            MetaObject mo(tuple.at(2)), mop(tuple.at(3));
            Object o(tuple.at(4)), op(tuple.at(5));
            Binary b(tuple.at(6));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));
            assert((0 <= mop) && (mop < num_meta_objects_));
            assert((0 <= b) && (b < num_static_binary_predicates_));

            if( (mo < mop) &&
                (o < objects_per_layer_.at(tr_layer_.at(t))) &&
                (op < objects_per_layer_.at(tr_layer_.at(t))) ) {
                Layer layer = tr_layer_.at(t);
                assert((0 <= layer) && (layer < num_layers_));
                assert((transitions_per_layer_.at(layer).front() <= t) && (t <= transitions_per_layer_.at(layer).back()));

                // map(t,a) & mapt(t,mo,o) & mapt(t,mop,op) & binary(b,a,mo,mop) => s(t.layer,b,o,op)
                add_implication({ 1 + map(t, a), 1 + mapt(t, mo, o), 1 + mapt(t, mop, op), 1 + binary(b, a, mo, mop) }, { 1 + s(layer, b, o, op) });
            }
        };

        SAT::VarSet tmp(p_transitions_, p_actions_, p_meta_objects_, p_meta_objects_, p_objects_, p_objects_, p_binary_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (Symmetries) Ordered arities for atoms:
    // arity(p,i) => OR { arity(p-1,j) : 0 <= j <= i }
    void build_formulas_arity_then_arity(std::ostream &os) {
        start_group("arity(p,i) => OR { arity(pp,j) : 1 + pp = p, 0 <= j <= i }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Atom p(tuple.at(0)), pp(p - 1);
            Arity i(tuple.at(1));
            assert((0 <= p) && (p < num_atoms_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (0 <= pp) && (0 <= i) ) {
                assert((0 <= pp) && (pp < num_atoms_));

                // arity(p,i) => OR { arity(pp,j) : 1 + pp = p, 0 <= j <= i }
                SAT::Implication IP({ 1 + arity(p, i) }, { });
                for( Arity j = Arity(0); j <= i; ++j )
                    IP.add_consequent(1 + arity(pp, j));
                add_implication(IP);
            }
        };

        arity.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (Symmetries) Ordered use of action arguments:
    // args(a,mo) => args(a,mop) [ 1 + mop = mo ]
    void build_formulas_args_then_args(std::ostream &os) {
        start_group("args(a,mo) => args(a,mop) [ 1 + mop = mo ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Action a(tuple.at(0));
            MetaObject mo(tuple.at(1)), mop(mo - 1);
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= mo) && (mo < num_meta_objects_));

            if( mop >= 0 ) {
                assert((0 <= mop) && (mop < num_meta_objects_));

                // args(a,mo) => args(a,mop) [ 1 + mop = mo ]
                add_implication({ 1 + args(a, mo) }, { 1 + args(a, mop) });
            }
        };

        SAT::VarSet tmp(p_actions_, p_meta_objects_);
        tmp.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (Symmetries) Definition of ord(o,k,i,s):
    // ord(o,k,i,s) <=> ground(k,i,o) & phi(k,s)
    void build_formulas_ord_iff_ground2_phi(std::ostream &os) {
        start_group("ord(o,k,i,s) <=> ground(k,i,o) & phi(k,s)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Object o(tuple.at(0));
            Feature k(tuple.at(1));
            Arity i(tuple.at(2));
            State s(tuple.at(3));
            assert((0 <= k) && (k < num_features_));
            assert((0 <= i) && (i < 1 + max_arity_));

            if( (o < objects_per_layer_.at(f_layer_.at(k))) && (0 < i) && (f_layer_.at(k) == s_layer_.at(s)) ) {
                Layer layer = f_layer_.at(k);
                assert((0 <= layer) && (layer < num_layers_));
                assert((states_per_layer_.at(layer).front() <= s) && (s <= states_per_layer_.at(layer).back()));

                // ord(o,k,i,s) => ground(k,i,o)
                // ord(o,k,i,s) => phi(k,s)
                // ord(o,k,i,s) <= ground(k,i,o) & phi(k,s)
                add_implication({ 1 + ord(o, k, i, s) }, { 1 + ground2(k, i, o) });
                add_implication({ 1 + ord(o, k, i, s) }, { 1 + phi(k, s) });
                add_implication({ 1 + ground2(k, i, o), 1 + phi(k, s) }, { 1 + ord(o, k, i, s) });
            }
        };

        ord.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // (Symmetries) Ordered objects at each layer
    // <non-strict-lexicographic-ordering-objects-in-layers>
    void build_formulas_lex_ordering_objects_in_layers(std::ostream &os) {
        start_group("<non-strict-lexicographic-ordering-objects-in-layers>");

        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            std::string prefix(std::string("objects-in-layer-") + std::to_string(layer));
            std::vector<std::vector<int> > lvectors;

            for( Object o = Object(0); o < objects_per_layer_.at(layer); ++o ) {
                std::vector<int> ords;
                for( int j1 = 0; j1 < int(states_per_layer_.at(layer).size()); ++j1 ) {
                    State s = states_per_layer_.at(layer).at(j1);
                    for( int j2 = 0; j2 < int(features_per_layer_.at(layer).size()); ++j2 ) {
                        Feature k = Feature(features_per_layer_.at(layer).at(j2));
                        for( Arity i = Arity(1); i < 1 + max_arity_; ++i )
                            ords.push_back(-(1 + ord(o, k, i, s)));
                    }
                }
                lvectors.emplace_back(std::move(ords));
            }

            if( !lvectors.empty() ) lex_ordering(prefix, lvectors, symmetries("strict-ordering-objects-in-layers"));
        }

        end_group(os);
    }

    // FULL (default) ENCODING OF GROUND ACTIONS:

    // Explanation of non-existing grounded actions:
    // -gtuple(l,a,<tuple>) => OR { -args(a,mo_i) : o_i > 0 } v OR { U(l,u,a,mo_i,o_i) : u, i } v OR { B(l,b,a,mo_i,mo_j,o_i,o_j) : b, i < j }
    void build_formulas_gtuple_then_args_U_B(std::ostream &os) {
        start_group("(full) -gtuple(l,a,<tuple>) => OR { -args(a,mo_i) : o_i > 0 } v OR { U(l,u,a,mo_i,o_i) : u, i } v OR { B(l,b,a,mo_i,mo_j,o_i,o_j) : b, i < j }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Layer layer(tuple.at(0));

            // enforce constraint only when handling a "complete" layer
            if( random_transitions_per_layer_.at(layer) == -1 ) {
                assert(int(tuple.size()) == 2 + num_meta_objects_);
                Action a = Action(tuple.at(1));
                assert((0 <= layer) && (layer < num_layers_));
                assert((0 <= a) && (a < num_actions_));

                for( int i = 0; i < num_meta_objects_; ++i ) {
                    assert(0 <= tuple.at(2 + i));
                    if( tuple.at(2 + i) >= objects_per_layer_.at(layer) )
                        return;
                }

                // -gtuple(l,a,<tuple>) => OR { -args(a,mo_i) : o_i > 0 } v OR { U(l,u,a,mo_i,o_i) : u, i } v OR { B(l,b,a,mo_i,mo_j,o_i,o_j) : b, i < j }
                SAT::Implication IP({ -(1 + gtuple(tuple)) }, { });

                // OR { -args(a,mo_i) : o_i > 0 }
                for( int i = 0; i < num_meta_objects_; ++i ) {
                    MetaObject mo_i = MetaObject(i);
                    if( tuple.at(2 + i) > 0 )
                        IP.add_consequent(-(1 + args(a, mo_i)));
                }

                // OR { U(l,u,a,mo_i,o_i) : u, i }
                for( int i = 0; i < num_meta_objects_; ++i ) {
                    MetaObject mo_i = MetaObject(i);
                    Object o_i = Object(tuple.at(2 + i));
                    for( Unary u = Unary(0); u < num_static_unary_predicates_; ++u )
                        IP.add_consequent(1 + U(layer, u, a, mo_i, o_i));
                }

                // OR { B(l,b,a,mo_i,mo_j,o_i,o_j) : b, i < j }
                for( int i = 0; i < num_meta_objects_; ++i ) {
                    MetaObject mo_i = MetaObject(i);
                    Object o_i = Object(tuple.at(2 + i));
                    for( int j = 1 + i; j < num_meta_objects_; ++j ) {
                        MetaObject mo_j = MetaObject(j);
                        Object o_j = Object(tuple.at(2 + j));
                        for( Binary b = Binary(0); b < num_static_binary_predicates_; ++b )
                            IP.add_consequent(1 + B(layer, b, a, mo_i, mo_j, o_i, o_j));
                    }
                }

                add_implication(IP);
            }
        };

        gtuple.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of gtuple(l,a,<tuple>) and G(t,a,<tuple>):
    // gtuple(l,a,<tuple>) => OR { G(t,a,<tuple>) : t.layer = l }
    void build_formulas_gtuple_then_G(std::ostream &os) {
        start_group("(full) gtuple(l,a,<tuple>) => OR { G(t,a,<tuple>) : t.layer = l }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Layer layer = Layer(tuple.at(0));
            Action a = Action(tuple.at(1));
            assert((0 <= layer) && (layer < num_layers_));
            assert((0 <= a) && (a < num_actions_));

            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(2 + i));
                if( tuple.at(2 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // gtuple(l,a,<tuple>) => OR { G(t,a,<tuple>) : t.layer = l }
            SAT::Implication IP({ 1 + gtuple(tuple) }, { });
            std::vector<int> args(tuple);
            for( int i = 0; i < int(transitions_per_layer_.at(layer).size()); ++i ) {
                Transition t = Transition(transitions_per_layer_.at(layer).at(i));
                assert(tr_layer_.at(t) == layer);
                args[0] = t;
                IP.add_consequent(1 + G(args));
            }
            add_implication(IP);
        };

        gtuple.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of gtuple(l,a,<tuple>) and G(t,a,<tuple>):
    // G(t,a,<tuple>) => gtuple(t.layer,a,<tuple>)
    void build_formulas_G_then_gtuple(std::ostream &os) {
        start_group("(full) G(t,a,<tuple>) => gtuple(t.layer,a,<tuple>)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Transition t = Transition(tuple.at(0));
            Action a = Action(tuple.at(1));
            assert((0 <= t) && (t < num_transitions_));
            assert((0 <= a) && (a < num_actions_));

            Layer layer = tr_layer_.at(t);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(2 + i));
                if( tuple.at(2 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // G(t,a,<tuple>) => gtuple(t.layer,a,<tuple>)
            std::vector<int> args(tuple);
            args[0] = layer;
            add_implication({ 1 + G(tuple) }, { 1 + gtuple(args) });
        };

        G.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Definition of gtuple(l,a,<tuple>) and G(t,a,<tuple>):
    // G(t,a,<tuple>) => map(t,a) & AND { mapt(t,mo_i,o_i) : o_i > 0 } AND { args(a,mo_i) => mapt(t,mo_i,o_i) : o_i = 0 }
    void build_formulas_G_then_map_mapt_args_then_mapt(std::ostream &os) {
        start_group("(full) G(t,a,<tuple>) => map(t,a) & AND { mapt(t,mo_i,o_i) : o_i > 0 } & AND { args(a,mo_i) => mapt(t,mo_i,o_i) : o_i = 0 }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Transition t = Transition(tuple.at(0));
            Action a = Action(tuple.at(1));
            Layer layer = tr_layer_.at(t);
            assert((0 <= t) && (t < num_transitions_));
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= layer) && (layer < num_layers_));

            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(2 + i));
                if( tuple.at(2 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // G(t,a,<tuple>) => map(t,a)
            add_implication({ 1 + G(tuple) }, { 1 + map(t, a) });

            // G(t,a,<tuple>) => AND { mapt(t,mo_i,o_i) : o_i > 0 }
            for( int i = 0; i < num_meta_objects_; ++i ) {
                MetaObject mo_i = MetaObject(i);
                Object o_i = Object(tuple.at(2 + i));
                if( o_i > 0 )
                    add_implication({ 1 + G(tuple) }, { 1 + mapt(t, mo_i, o_i) });
            }

            // G(t,a,<tuple>) => AND { args(a,mo_i) => mapt(t,mo_i,o_i) : o_i = 0 }
            for( int i = 0; i < num_meta_objects_; ++i ) {
                MetaObject mo_i = MetaObject(i);
                Object o_i = Object(tuple.at(2 + i));
                if( o_i == 0 )
                    add_implication({ 1 + G(tuple), 1 + args(a, mo_i) }, { 1 + mapt(t, mo_i, o_i) });
            }
        };

        G.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Explanation of existing grounded actions:
    // AT-MOST-1 { G(t,a,<tuple>) : t.src = s } [ s, a(<tuple>) ]
    void build_formulas_at_most_1_G(std::ostream &os) {
        start_group("(full) AT-MOST-1 { G(t,a,<tuple>) : t.src = s } [ s, a(<tuple>) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            Layer layer(tuple.at(0));
            std::vector<int> ntuple(tuple);
            for( int i = 0; i < int(states_per_layer_.at(layer).size()); ++i ) {
                State s = states_per_layer_.at(layer)[i];
                std::vector<int> literals;
                for( int j = 0; j < int(transitions_per_layer_.at(layer).size()); ++j ) {
                    ntuple[0] = transitions_per_layer_.at(layer)[j];
                    if( tr_src_.at(ntuple[0]) == s ) literals.push_back(1 + G(ntuple));
                }
                assert(!literals.empty());
                at_most_1(std::string("G-at-most-1(") + std::to_string(s) + ")", literals);
            }
        };

        gtuple.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Explanation of existing grounded actions:
    // EXACTLY-1 { G(t,a,<tuple>) : a, <tuple> }
    void build_formulas_exactly_1_G(std::ostream &os) {
        start_group("(full) EXACTLY-1 { G(t,a,<tuple>) : a, <tuple> }");

        for( Transition t = Transition(0); t < num_transitions_; ++t ) {
            std::vector<int> literals;
            auto foo = [this, t, &literals](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                if( tuple.at(0) == t ) literals.push_back(1 + G(tuple));
            };
            G.enumerate_vars_from_multipliers(foo);
            exactly_1(std::string("G-exactly-1(") + std::to_string(t) + ")", literals);
        }

        end_group(os);
    }

    // Explanation of existing grounded actions:
    // AT-LEAST-1 { G(t,a,<tuple>) : a, <tuple> }
    void build_formulas_at_least_1_G(std::ostream &os) {
        start_group("(full) AT-LEAST-1 { G(t,a,<tuple>) : a, <tuple> }");

        for( Transition t = Transition(0); t < num_transitions_; ++t ) {
            std::vector<int> literals;
            auto foo = [t, &literals, this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                if( tuple.at(0) == t ) literals.push_back(1 + G(tuple));
            };
            G.enumerate_vars_from_multipliers(foo);
            assert(!literals.empty());
            at_least_1(std::string("G-at-least-1(") + std::to_string(t) + ")", literals);
        }

        end_group(os);
    }

    // Applicable actions must be applied:
    // G(t,a,<tuple>) => appl(a,<tuple>,t.src)
    void build_formulas_G_then_appl2(std::ostream &os) {
        start_group("(full) G(t,a,<tuple>) => appl(a,<tuple>,t.src)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Transition t = Transition(tuple.front());
            Action a = Action(tuple.at(1));
            assert((0 <= t) && (t < num_transitions_));
            assert((0 <= a) && (a < num_actions_));

            Layer layer = tr_layer_.at(t);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(2 + i));
                if( tuple.at(2 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // make new tuple for appl(a,<tuple>,t.src)
            std::vector<int> ntuple{ a };
            ntuple.insert(ntuple.end(), &tuple[2], &tuple[2 + num_meta_objects_]);
            ntuple.push_back(tr_src_.at(t));

            // G(t,a,<tuple>) => appl(a,<tuple>,t.src)
            add_implication({ 1 + G(tuple) }, { 1 + appl2(ntuple) });
        };

        G.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // appl(a,<tuple>,s) => OR { G(t,a,<tuple>) : t.src = s }
    void build_formulas_appl2_then_G(std::ostream &os) {
        start_group("(full) appl(a,<tuple>,s) => OR { G(t,a,<tuple>) : t.src = s }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Action a = Action(tuple.front());
            State s = State(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= s) && (s < num_states_));

            Layer layer = s_layer_.at(s);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(1 + i));
                if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // appl(a,<tuple>,s) => OR { G(t,a,<tuple>) : t.src = s }
            SAT::Implication IP({ 1 + appl2(tuple) }, { });
            for( int t = 0; t < num_transitions_; ++t ) {
                if( tr_src_.at(t) == s ) {
                    assert(tr_layer_.at(t) == layer);
                    std::vector<int> ntuple{ t };
                    ntuple.insert(ntuple.end(), &tuple[0], &tuple[1 + num_meta_objects_]);
                    IP.add_consequent(1 + G(ntuple));
                }
            }
            add_implication(IP);
        };

        appl2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // -appl(a,<tuple>,s) => -gtuple(s.layer,a,<tuple>) v OR { violated0(a,<tuple>,s,k) : k } v OR { violated1(a,<tuple>,s,k) : k }
    void build_formulas_appl2_then_gtuple_violated0_violated1(std::ostream &os) {
        start_group("(full) -appl(a,<tuple>,s) => -gtuple(s.layer,a,<tuple>) v OR { violated0(a,<tuple>,s,k) : k } v OR { violated1(a,<tuple>,s,k) : k }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            Action a = Action(tuple.front());
            State s = State(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= s) && (s < num_states_));

            Layer layer = s_layer_.at(s);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(1 + i));
                if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // -appl(a,<tuple>,s) => -gtuple(s.layer,a,<tuple>) v OR { violated0(a,<tuple>,s,k) : k } v OR { violated1(a,<tuple>,s,k) : k }
            SAT::Implication IP({ -(1 + appl2(tuple)) }, { });

            // -gtuple(s.layer,a,<tuple>)
            std::vector<int> ntuple{ layer, a };
            ntuple.insert(ntuple.end(), &tuple[1], &tuple[1 + num_meta_objects_]);
            IP.add_consequent(-(1 + gtuple(ntuple)));

            // OR { violated0(a,<tuple>,s,k) : k }
            for( Feature k = Feature(0); k < num_features_; ++k ) {
                if( f_layer_.at(k) == layer ) {
                    std::vector<int> ntuple(tuple);
                    ntuple.push_back(k);
                    IP.add_consequent(1 + violated0(ntuple));
                }
            }

            // OR { violated1(a,<tuple>,s,k) : k }
            for( Feature k = Feature(0); k < num_features_; ++k ) {
                if( f_layer_.at(k) == layer ) {
                    std::vector<int> ntuple(tuple);
                    ntuple.push_back(k);
                    IP.add_consequent(1 + violated1(ntuple));
                }
            }

            add_implication(IP);
        };

        appl2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // violated0(a,<tuple>,s,k) => phi(k,s) & OR { pre0eq(a,<tuple>,k,mk) : mk }
    void build_formulas_violated0_then_phi_pre0eq(const std::vector<int> &tuple) {
    }
    void build_formulas_violated0_then_phi_pre0eq(std::ostream &os) {
        start_group("(full) violated0(a,<tuple>,s,k) => phi(k,s) & OR { pre0eq(a,<tuple>,k,mk) : mk }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 3 + num_meta_objects_);
            Action a = Action(tuple.front());
            State s = State(tuple.at(1 + num_meta_objects_));
            Feature k = Feature(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= s) && (s < num_states_));
            assert((0 <= k) && (k < num_features_));

            if( s_layer_.at(s) == f_layer_.at(k) ) {
                Layer layer = s_layer_.at(s);
                for( int i = 0; i < num_meta_objects_; ++i ) {
                    assert(0 <= tuple.at(1 + i));
                    if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                        return;
                }

                // violated0(a,<tuple>,s,k) => phi(k,s)
                add_implication({ 1 + violated0(tuple) }, { 1 + phi(k, s) });

                // violated0(a,<tuple>,s,k) => OR { pre0eq(a,<tuple>,k,mk) : mk }
                SAT::Implication IP({ 1 + violated0(tuple) }, { });
                std::vector<int> ntuple(tuple);
                ntuple[1 + num_meta_objects_] = k;
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
                    ntuple[2 + num_meta_objects_] = mk;
                    IP.add_consequent(1 + pre0eq(ntuple));
                }
                add_implication(IP);
            }
        };

        violated0.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // violated1(a,<tuple>,s,k) => -phi(k,s) & OR { pre1eq(a,<tuple>,k,mk) : mk }
    void build_formulas_violated1_then_phi_pre1eq(std::ostream &os) {
        start_group("(full) violated1(a,<tuple>,s,k) => -phi(k,s) & OR { pre1eq(a,<tuple>,k,mk) : mk }");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 3 + num_meta_objects_);
            Action a = Action(tuple.front());
            State s = State(tuple.at(1 + num_meta_objects_));
            Feature k = Feature(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= s) && (s < num_states_));
            assert((0 <= k) && (k < num_features_));

            if( s_layer_.at(s) == f_layer_.at(k) ) {
                Layer layer = s_layer_.at(s);
                for( int i = 0; i < num_meta_objects_; ++i ) {
                    assert(0 <= tuple.at(1 + i));
                    if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                        return;
                }

                // violated1(a,<tuple>,s,k) => -phi(k,s)
                add_implication({ 1 + violated1(tuple) }, { -(1 + phi(k, s)) });

                // violated1(a,<tuple>,s,k) => OR { pre1eq(a,<tuple>,k,mk) : mk }
                SAT::Implication IP({ 1 + violated1(tuple) }, { });
                std::vector<int> ntuple(tuple);
                ntuple[1 + num_meta_objects_] = k;
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
                    ntuple[2 + num_meta_objects_] = mk;
                    IP.add_consequent(1 + pre1eq(ntuple));
                }
                add_implication(IP);
            }
        };

        violated1.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // pre0eq(a,<tuple>,k,mk) => pre0(a,mk) & eq(<tuple>,mk,k)
    void build_formulas_pre0eq_then_pre0_eq(std::ostream &os) {
        start_group("(full) pre0eq(a,<tuple>,k,mk) => pre0(a,mk) & eq(<tuple>,mk,k)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 3 + num_meta_objects_);

            Action a = Action(tuple.front());
            Feature k = Feature(tuple.at(1 + num_meta_objects_));
            MetaFeature mk = MetaFeature(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= k) && (k < num_features_));
            assert((0 <= mk) && (mk < num_meta_features_));

            Layer layer = f_layer_.at(k);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(1 + i));
                if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // pre0eq(a,<tuple>,k,mk) => pre0(a,mk)
            add_implication({ 1 + pre0eq(tuple) }, { 1 + pre0(a, mk) });

            // pre0eq(a,<tuple>,k,mk) => eq(<tuple>,mk,k)
            std::vector<int> ntuple(&tuple[1], &tuple[1 + num_meta_objects_]);
            ntuple.push_back(mk);
            ntuple.push_back(k);
            assert(int(ntuple.size()) == 2 + num_meta_objects_);
            add_implication({ 1 + pre0eq(tuple) }, { 1 + eq2(ntuple) });
        };

        pre0eq.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // pre1eq(a,<tuple>,k,mk) => pre1(a,mk) & eq(<tuple>,mk,k)
    void build_formulas_pre1eq_then_pre1_eq(const std::vector<int> &tuple) {
    }
    void build_formulas_pre1eq_then_pre1_eq(std::ostream &os) {
        start_group("(full) pre1eq(a,<tuple>,k,mk) => pre1(a,mk) & eq(<tuple>,mk,k)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 3 + num_meta_objects_);

            Action a = Action(tuple.front());
            Feature k = Feature(tuple.at(1 + num_meta_objects_));
            MetaFeature mk = MetaFeature(tuple.back());
            assert((0 <= a) && (a < num_actions_));
            assert((0 <= k) && (k < num_features_));
            assert((0 <= mk) && (mk < num_meta_features_));

            Layer layer = f_layer_.at(k);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(1 + i));
                if( tuple.at(1 + i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // pre1eq(a,<tuple>,k,mk) => pre1(a,mk)
            add_implication({ 1 + pre1eq(tuple) }, { 1 + pre1(a, mk) });

            // pre1eq(a,<tuple>,k,mk) => eq(<tuple>,mk,k)
            std::vector<int> ntuple(&tuple[1], &tuple[1 + num_meta_objects_]);
            ntuple.push_back(mk);
            ntuple.push_back(k);
            assert(int(ntuple.size()) == 2 + num_meta_objects_);
            add_implication({ 1 + pre1eq(tuple) }, { 1 + eq2(ntuple) });
        };

        pre1eq.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // eq(<tuple>,mk,k) => [ atom(mk,p) <=> ground(k,p) ]
    void build_formulas_eq_atom1_then_ground1(std::ostream &os) {
        start_group("(full) eq(<tuple>,mk,k) => [ atom(mk,p) <=> ground(k,p) ]");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            MetaFeature mk = MetaFeature(tuple.at(num_meta_objects_));
            Feature k = Feature(tuple.back());
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= k) && (k < num_features_));

            Layer layer = f_layer_.at(k);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(i));
                if( tuple.at(i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // eq(<tuple>,mk,k) & atom(mk,p) => ground(k,p)
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                add_implication({ 1 + eq2(tuple), 1 + atom1(mk, p) }, { 1 + ground1(k, p) });

            // eq(<tuple>,mk,k) & ground(k,p) => atom(mk,p)
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                add_implication({ 1 + eq2(tuple), 1 + ground1(k, p) }, { 1 + atom1(mk, p) });
        };

        eq2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Applicable actions must be applied:
    // eq(<tuple>,mk,k) & atom(mk,i,mo_j) => ground(k,i,o_j)
    void build_formulas_eq_atom2_then_ground2(std::ostream &os) {
        start_group("(full) eq(<tuple>,mk,k) & atom(mk,i,mo_j) => ground(k,i,o_j)");

        auto foo = [this](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            assert(int(tuple.size()) == 2 + num_meta_objects_);

            MetaFeature mk = MetaFeature(tuple.at(num_meta_objects_));
            Feature k = Feature(tuple.back());
            assert((0 <= mk) && (mk < num_meta_features_));
            assert((0 <= k) && (k < num_features_));

            Layer layer = f_layer_.at(k);
            for( int i = 0; i < num_meta_objects_; ++i ) {
                assert(0 <= tuple.at(i));
                if( tuple.at(i) >= objects_per_layer_.at(layer) )
                    return;
            }

            // eq(<tuple>,mk,k) & atom(mk,i,mo_j) => ground(k,i,o_j)
            for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    add_implication({ 1 + eq2(tuple), 1 + atom2(mk, i, mo) }, { 1 + ground2(k, i, tuple.at(mo)) });
            }
        };

        eq2.enumerate_vars_from_multipliers(foo);
        end_group(os);
    }

    // Units for verification of meta-layer:
    void build_formulas_verification_meta_layer(std::ostream &os) {
        start_group("units for verification of meta-layer");

        for( Layer layer = Layer(0); layer < num_layers_; ++layer ) {
            Feature k = features_per_layer_.at(layer).front();
            for( Atom p = Atom(0); p < num_atoms_; ++p ) {
                SAT::VarSet tmp;
                std::vector<int> objects(objects_per_layer_.at(layer));
                std::iota(objects.begin(), objects.end(), 0);
                for( Arity i = Arity(0); i < arity_for_atoms_.at(p); ++i )
                    tmp.fill_multipliers(objects);

                auto foo = [this, layer, &k, p](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                    assert(arity_for_atoms_.at(p) == int(tuple.size()));
                    assert(k <= features_per_layer_.at(layer).back());

                    // ground(k,p)
                    add_unit(1 + ground1(k, p));

                    // ground(k,i,o)
                    for( Arity i = Arity(1); i < 1 + arity_for_atoms_.at(p); ++i ) {
                        Object o(tuple.at(i - 1));
                        add_unit(1 + ground2(k, i, o));
                    }

                    // increment feature
                    ++k;
                };
                tmp.enumerate_vars_from_multipliers(foo);
            }
        }

        end_group(os);
    }

    void build_meta_layer(std::ostream &os) {
        os << green("Basic definitions:", !options_.disable_colors_) << std::endl;
        build_formulas_using1_iff_using2(os);
        build_formulas_using2_iff_pre0_pre1_eff0_eff1(os);

        os << green("Consistent preconditions, effects, and labeling:", !options_.disable_colors_) << std::endl;
        build_formulas_then_pre0_pre1(os);
        build_formulas_then_eff0_eff1(os);
        build_formulas_at_most_1_label(os);

        os << green("Effects are non-redundant:", !options_.disable_colors_) << std::endl;
        build_formulas_eff_then_pre(os);

        if( regularizer("disjoint-meta-features") ) {
            os << yellow("(Regularizer)", !options_.disable_colors_) << green(" Disjoint meta-features in actions:", !options_.disable_colors_) << std::endl;
            build_formulas_at_most_1_using2(os);
        }

        os << green("Consistent arities for atoms:", !options_.disable_colors_) << std::endl;
        build_formulas_exactly_1_arity(os);

        os << green("Mapping of meta-features to atoms:", !options_.disable_colors_) << std::endl;
        build_formulas_exactly_1_atom1(os);
        build_formulas_at_most_1_atom2(os);
        build_formulas_atom1_atom2_then_arity(os);
        build_formulas_atom1_arity_then_atom2(os);
        build_formulas_atom1_arity_then_atom2_2(os);
        build_formulas_atom2(os);

        os << green("1-1 mapping of meta-features to atoms:", !options_.disable_colors_) << std::endl;
        build_formulas_lex_ordering_meta_layer(os);

        os << green("Atoms must be non-static:", !options_.disable_colors_) << std::endl;
        build_formulas_non_static0_non_static1(os);
        build_formulas_non_static0_then_atom_pre1_eff0(os);
        build_formulas_non_static1_then_atom_pre0_eff1(os);

        os << green("Relevant arguments for actions:", !options_.disable_colors_) << std::endl;
        build_formulas_using2_atom2_then_args(os);
        build_formulas_args_then_relevant(os);
        build_formulas_relevant_iff_using2_atom2(os);

        os << green("Arity for actions and atoms:", !options_.disable_colors_) << std::endl;
        build_formulas_args(os);
        build_formulas_arity(os);

        if( regularizer("exact-arities") ) {
            os << yellow("(Regularizer)", !options_.disable_colors_) << green(" Exact arities:", !options_.disable_colors_) << std::endl;
            build_formulas_exact_args(os);
            build_formulas_exact_arity(os);
        }

        os << green("Static predicates on relevant arguments:", !options_.disable_colors_) << std::endl;
        build_formulas_unary_then_args(os);
        build_formulas_binary_then_args_args(os);
    }

    void build_layer(std::ostream &os) {
        os << green("Consistent mapping of transitions to actions:", !options_.disable_colors_) << std::endl;
        build_formulas_exactly_1_map(os);

        os << green("Consistent mapping of features to meta-features in actions:", !options_.disable_colors_) << std::endl;
        build_formulas_at_most_1_mapf_by_k(os);
        build_formulas_at_most_1_mapf_by_mk(os);

        os << green("Consistency between map, mapf, labeling, and using:", !options_.disable_colors_) << std::endl;
        build_formulas_map_then_label(os);
        build_formulas_map_mapf_then_using2(os);
        build_formulas_map_using2_then_mapf(os);

        os << green("Definition of free(k,t,a):", !options_.disable_colors_) << std::endl;
        build_formulas_map_mapf_then_free(os);
        build_formulas_map_mapf_then_eff0_eff1_iff_free(os);

        os << green("Transitions + inertia:", !options_.disable_colors_) << std::endl;
        build_formulas_transitions(os);
        build_formulas_transitions_inertia(os);

        if( encoding("applicable-actions-alt") ) {
            os << green("Applicable actions must be applied:", !options_.disable_colors_) << std::endl;
            build_formulas_appl_then_mapeq(os);
            build_formulas_mapeq_then_map_eq(os);
            build_formulas_eq_then_mapf_iff_mapf(os);
        }

        os << green("Definition of g(k,s,t):", !options_.disable_colors_) << std::endl;
        build_formulas_def_g(os);

        os << green("Separate different states using features:", !options_.disable_colors_) << std::endl;
        build_formulas_separate_states(os);

        if( encoding("applicable-actions-alt") ) {
            os << green("Definition of appl(a,t,s):", !options_.disable_colors_) << std::endl;
            build_formulas_map_Z0_Z1_then_appl(os);

            os << green("Definition of Z0(t,k,a,s):", !options_.disable_colors_) << std::endl;
            build_formulas_phi_then_Z0(os);
            build_formulas_X0_then_Z0(os);
            build_formulas_X0_then_pre0_mapf(os);

            os << green("Definition of Z1(t,k,a,s):", !options_.disable_colors_) << std::endl;
            build_formulas_phi_then_Z1(os);
            build_formulas_X1_then_Z1(os);
            build_formulas_X1_then_pre1_mapf(os);
        }

        os << green("Mapping of features to grounded atoms:", !options_.disable_colors_) << std::endl;
        build_formulas_exactly_1_ground1(os);
        build_formulas_at_most_1_ground2(os);
        build_formulas_ground1_ground2_then_arity(os);
        build_formulas_ground1_arity_then_ground2(os);
        build_formulas_ground1_arity_then_ground2_2(os);
        build_formulas_ground2(os);

        if( encoding("features-to-actions-map-alt1") ) {
            os << green("1-1 mapping of features to ground atoms:", !options_.disable_colors_) << std::endl;
            ENC1_build_formulas_ground1_ground1_arity_then_gdiff(os);
            ENC1_build_formulas_gdiff_then_ground2_ground2(os);
            ENC1_build_formulas_gdiff_then_ground2_ground2_2(os);
            ENC1_build_formulas_gdiff_then_gdiff(os);
        } else if( encoding("features-to-actions-map-alt2") ) {
            os << green("1-1 mapping of features to ground atoms:", !options_.disable_colors_) << std::endl;
            ENC2_build_formulas_arity_then_ground1_ground1(os);
            ENC2_build_formulas_gdiff(os);
            ENC2_build_formulas_gdiff_then_ground1_ground1(os);
            ENC2_build_formulas_gdiff_then_ground2_ground2(os);
            ENC2_build_formulas_gdiff_then_ground2_ground2_2(os);
            // [DISABLED] ENC2_build_formulas_at_most_1_gdiff(os);
        } else if( encoding("features-to-actions-map-alt3") ) {
            os << green("1-1 mapping of features to ground atoms:", !options_.disable_colors_) << std::endl;
            ENC3_build_formulas_gdiff1_iff_ground1_ground1(os);
            ENC3_build_formulas_gdiff1_arity_then_gdiff2(os);
            ENC3_build_formulas_gdiff2_then_gdiff1_arity(os);
            ENC3_build_formulas_gdiff2_then_ground2_ground2(os);
        }

        if( some_features_to_actions_map_encoding() ) {
            os << green("1-1 mapping of features to ground atoms:", !options_.disable_colors_) << std::endl;
            build_formulas_lex_ordering_features_in_layers(os);
        }

        os << green("Consistent mapping between features and meta-features:", !options_.disable_colors_) << std::endl;
        build_formulas_mapf_then_atom1_iff_ground1(os);
        build_formulas_mapf_atom2_then_ground2(os);
        build_formulas_mapf_ground2_then_atom2(os);

        os << green("Definition of U(l,u,a,mo,o) and B(l,b,a,mo,mop,o,op):", !options_.disable_colors_) << std::endl;
        build_formulas_U_iff_unary_r(os);
        build_formulas_B_iff_binary_s(os);

        os << green("Definition of mapt(t,mo,o):", !options_.disable_colors_) << std::endl;
        build_formulas_at_most_1_mapt(os);
        build_formulas_map_args_then_mapt(os);
        build_formulas_map_mapt_then_args(os);

        os << green("Cross consistency between schemas and transitions:", !options_.disable_colors_) << std::endl;
        //CHECK build_formulas_mapf_atom2_then_ground2_iff_mapt(os);
        build_formulas_mapf_atom2_then_W(os);
        build_formulas_W_then_ground2_iff_mapt(os);

        os << green("One-hot encoding of variables:", !options_.disable_colors_) << std::endl;
        build_formulas_feature_sum(os);
        if( some_features_to_actions_map_encoding() && !encoding("features-to-actions-map-strict-lex-ordering") ) {
            build_formulas_negative_units_features(os);
            build_formulas_positive_units_features(os);
        }

        os << green("Explanation of existing grounded actions (common to different encodings):", !options_.disable_colors_) << std::endl;
        build_formulas_map_mapt_unary_then_r(os);
        build_formulas_map_mapt_mapt_binary_then_s(os);

        if( symmetries("ordering-arities-for-atoms") ) {
            os << yellow("(Symmetries)", !options_.disable_colors_) << green(" Ordered arities for atoms:", !options_.disable_colors_) << std::endl;
            build_formulas_arity_then_arity(os);
        }

        if( symmetries("ordering-action-arguments") ) {
            os << yellow("(Symmetries)", !options_.disable_colors_) << green(" Ordered use of action arguments:", !options_.disable_colors_) << std::endl;
            build_formulas_args_then_args(os);
        }

        if( ordering_objects_in_layers() && some_features_to_actions_map_encoding() ) {
            os << yellow("(Symmetries)", !options_.disable_colors_) << green(" Definition of ord(o,k,i,s):", !options_.disable_colors_) << std::endl;
            build_formulas_ord_iff_ground2_phi(os);

            os << yellow("(Symmetries)", !options_.disable_colors_) << green(" Ordered objects at each layer:", !options_.disable_colors_) << std::endl;
            build_formulas_lex_ordering_objects_in_layers(os);
        }

        // FULL (default) ENCODING OF GROUND ACTIONS:

        os << green("Explanation of non-existing grounded actions:", !options_.disable_colors_) << std::endl;
        build_formulas_gtuple_then_args_U_B(os);

        os << green("Definition of gtuple(l,a,<tuple>) and G(t,a,<tuple>):", !options_.disable_colors_) << std::endl;
        build_formulas_gtuple_then_G(os);
        build_formulas_G_then_gtuple(os);
        build_formulas_G_then_map_mapt_args_then_mapt(os);

        os << green("Explanation of existing grounded actions:", !options_.disable_colors_) << std::endl;
        build_formulas_at_most_1_G(os);
        build_formulas_exactly_1_G(os);
        // [DISABLED] build_formulas_at_least_1_G(os);
        //build_formulas_gtuple_unary_then_r(os);   // [SUBSUMED]
        //build_formulas_gtuple_binary_then_s(os);  // [SUBSUMED]

        if( encoding("applicable-actions-tuples") ) {
            os << green("Applicable actions must be applied:", !options_.disable_colors_) << std::endl;
            build_formulas_G_then_appl2(os);
            build_formulas_appl2_then_G(os);
            build_formulas_appl2_then_gtuple_violated0_violated1(os);
            build_formulas_violated0_then_phi_pre0eq(os);
            build_formulas_violated1_then_phi_pre1eq(os);
            build_formulas_pre0eq_then_pre0_eq(os);
            build_formulas_pre1eq_then_pre1_eq(os);
            build_formulas_eq_atom1_then_ground1(os);
            build_formulas_eq_atom2_then_ground2(os);
        }

        if( options_.verify_meta_layer_ ) {
            os << green("Units for verification of meta-layer:", !options_.disable_colors_) << std::endl;
            build_formulas_verification_meta_layer(os);
        }
    }

    void build_base() override {
        build_meta_layer(std::cout);
        build_layer(std::cout);
    }

    void build_rest() override {
        if( !options_.partial_assignment_ && !options_.disable_symmetry_handling_ ) {
            // assign labels to first set of actions
            assert(int(num_labels_) <= int(num_actions_));
            for( Label l = Label(0); l < num_labels_; ++l ) {
                SAT::Implication *IP = new SAT::Implication;
                IP->add_consequent(1 + label(Action(l), l));
                add_implication(IP);
            }
            std::cout << "#units for action labels=" << num_labels_ << std::endl;
        }

        // limit number ground actions per layer
        int start = num_implications();
        for( Layer l = Layer(0); l < num_layers_; ++l ) {
            if( options_.max_number_ground_actions_.find(l) != options_.max_number_ground_actions_.end() ) {
                int n = options_.max_number_ground_actions_.at(l);
                std::vector<int> literals;
                auto foo = [this, l, &literals](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                    if( tuple.at(0) == l )
                        literals.push_back(1 + gtuple(tuple));
                };
                gtuple.enumerate_vars_from_multipliers(foo);

                if( !literals.empty() && (n < int(literals.size())) ) {
                    std::vector<int> z;
                    sorting_network(std::string("bound-num-ground-actions(layer=") + std::to_string(l) + ")", literals, z);
                    assert(z.size() == literals.size());

                    for( int j = n; j < int(literals.size()); ++j ) {
                        SAT::Implication *IP = new SAT::Implication;
                        IP->add_consequent(-(1 + z.at(j)));
                        add_implication(IP);
                    }
                }
            }
        }
        if( num_implications() - start > 0 )
            std::cout << "#implications for bound on ground actions=" << num_implications() - start << std::endl;
    }

    void build_soft_theory() override {
        if( options_.weighted_sat_ ) {
            int weight = 0;
            auto foo = [this, &weight](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                int var = varset(tuple);
                SAT::Implication *IP = new SAT::Implication;
                IP->add_consequent(-(1 + var));
                add_soft_implication(weight, IP);
            };

#if 0
            std::cout << yellow("Minimizing extension of static predicates", !options_.disable_colors_) << std::endl;
            int start = num_soft_implications();
            weight = 1;
            r.enumerate_vars_from_multipliers(foo);
            std::cout << "-r(l,u,o) (weight=1): #soft-implications=" << num_soft_implications() - start << std::endl;
            start = num_soft_implications();
            weight = 2;
            s.enumerate_vars_from_multipliers(foo);
            std::cout << "-s(l,b,o,op) (weight=2): #soft-implications=" << num_soft_implications() - start << std::endl;
#else
            std::cout << yellow("Minimizing number of grounded actions", !options_.disable_colors_) << std::endl;
            int start = num_soft_implications();
            weight = 1;
            gtuple.enumerate_vars_from_multipliers(foo);
            std::cout << "-gtuple(l,a,<tuple>) (weight=1): #soft-implications=" << num_soft_implications() - start << std::endl;
#endif
        }
    }

  public:
    // action to string
    std::string action_to_string(Action a) const {
        return std::string("a") + std::to_string(a);
    }

    // action label
    std::string action_label(Action a) const {
        std::string desc;
        for( Label l = Label(0); l < num_labels_; ++l ) {
            if( model_.at(label(a, l)) )
                desc += label_as_string_.at(l);
        }
        return desc;
    }

    // meta-feature to string
    std::string meta_feature_to_string(MetaFeature mk) const {
        return std::string("mk") + std::to_string(mk);
    }

    // feature to string
    std::string feature_to_string(Feature k) const {
        std::string str(std::string("f") + std::to_string(k));
        return str;
    }

    // meta-object to string
    std::string meta_object_to_string(MetaObject mo) const {
        return std::string("?obj") + std::to_string(mo);
    }

    // object to string
    std::string object_to_string(Object o) const {
        return std::string("o") + std::to_string(o);
    }

    // atom schema to string
    std::string atom_schema_desc(MetaFeature mk) const {
        std::string desc;
        bool atom_found = false;
        for( Atom p = Atom(0); p < num_atoms_; ++p ) {
            if( model_[atom1(mk, p)] ) {
                atom_found = true;
                desc += std::string("p") + std::to_string(p) + "(";

                Arity p_arity = Arity(-1);
                for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                    if( model_[arity(p, i)] ) {
                        p_arity = i;
                        break;
                    }
                }
                assert(p_arity != -1);

                for( Arity i = Arity(1); i <= p_arity; ++i ) {
                    bool meta_object_found = false;
                    for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
                        if( model_[atom2(mk, i, mo)] ) {
                            meta_object_found = true;
                            if( i > 1 ) desc += ",";
                            desc += meta_object_to_string(mo);
                        }
                    }
                    assert(meta_object_found);
                }
                desc += ")";
            }
        }
        assert(atom_found);
        return desc;
    }

    // ground atom to string
    std::string ground_atom_desc(Feature k) const {
        std::string desc;
        bool atom_found = false;
        for( Atom p = Atom(0); p < num_atoms_; ++p ) {
            if( model_[ground1(k, p)] ) {
                atom_found = true;
                desc += std::string("p") + std::to_string(p) + "(";

                Arity p_arity = Arity(-1);
                for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                    if( model_[arity(p, i)] ) {
                        p_arity = i;
                        break;
                    }
                }
                assert(p_arity != -1);

                for( Arity i = Arity(1); i <= p_arity; ++i ) {
                    bool object_found = false;
                    for( Object o = Object(0); o < objects_per_layer_.at(f_layer_.at(k)); ++o ) {
                        if( model_[ground2(k, i, o)] ) {
                            object_found = true;
                            if( i > 1 ) desc += ",";
                            desc += object_to_string(o);
                        }
                    }
                    assert(object_found);
                }
                desc += ")";
            }
        }
        assert(atom_found);
        return desc;
    }

    // unary schema to string
    std::string unary_predicate_to_string(Unary u) const {
        return std::string("u") + std::to_string(u);
    }
    std::string unary_schema_to_string(Unary u, MetaObject mo) const {
        return unary_predicate_to_string(u) + "(" + meta_object_to_string(mo) + ")";
    }

    // binary schema to string
    std::string binary_predicate_to_string(Binary b) const {
        return std::string("b") + std::to_string(b);
    }
    std::string binary_schema_to_string(Binary b, MetaObject mo, MetaObject mop) const {
        return binary_predicate_to_string(b) + "(" + meta_object_to_string(mo) + "," + meta_object_to_string(mop) + ")";
    }

    // description of meta-feature
    std::string meta_feature_desc(MetaFeature mk) const {
        std::string desc(std::string("meta-feature ") + meta_feature_to_string(mk) + ":");
        desc += std::string(" atom=") + atom_schema_desc(mk);
        desc += ", used-by=";
        bool need_comma = false;
        for( Action a = Action(0); a < num_actions_; ++a ) {
            if( model_[using2(a, mk)] ) {
                desc += !need_comma ? "{" : ",";
                desc += action_to_string(a);
                need_comma = true;
            }
        }
        desc += need_comma ? "}" : "<unused>";
        return desc;
    }

    // constant feature?
    bool non_constant_feature(Feature k) const {
        Layer l = f_layer_.at(k);
        bool value = model_.at(phi(k, states_per_layer_[l].front()));
        for( int i = 0; i < int(states_per_layer_[l].size()); ++i ) {
            if( value != model_.at(phi(k, states_per_layer_[l].at(i))) )
                return true;
        }
        return false;
    }

    // description of meta-feature
    std::string feature_desc(Feature k) const {
        std::string desc(std::string("feature ") + feature_to_string(k) + ":");
        desc += std::string(" atom=") + ground_atom_desc(k);
        desc += ", extension=";
        for( int i = 0; i < int(states_per_layer_[f_layer_.at(k)].size()); ++i )
            desc += std::to_string(model_.at(phi(k, states_per_layer_[f_layer_.at(k)].at(i))));
        return desc;
    }

    // description of action in terms of meta-features
    std::string action_mk_desc(Action a) const {
        std::string desc(std::string("action ") + action_to_string(a) + ":");
        desc += std::string(" label=") + action_label(a);

        // precondition
        desc += ", pre={";
        bool need_comma = false;
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            if( model_.at(pre0(a, mk)) ) {
                if( need_comma ) desc += ",";
                desc += std::string("-") + meta_feature_to_string(mk);
                need_comma = true;
            } else if( model_.at(pre1(a, mk)) ) {
                if( need_comma ) desc += ",";
                desc += meta_feature_to_string(mk);
                need_comma = true;
            }
        }
        desc += "}";

        // effect
        desc += ", eff={";
        need_comma = false;
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            if( model_.at(eff0(a, mk)) ) {
                if( need_comma ) desc += ",";
                desc += std::string("-") + meta_feature_to_string(mk);
                need_comma = true;
            }
            if( model_.at(eff1(a, mk)) ) {
                if( need_comma ) desc += ",";
                desc += meta_feature_to_string(mk);
                need_comma = true;
            }
        }
        desc += "}";

        // affected meta-features
        desc += ", aff={";
        need_comma = false;
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            if( model_.at(eff0(a, mk)) || model_.at(eff1(a, mk)) ) {
                if( need_comma ) desc += ",";
                desc += meta_feature_to_string(mk);
                need_comma = true;
            }
        }
        desc += "}";

        return desc;
    }

    // description of action schema
    std::string action_schema_desc(Action a) const {
        std::string desc(std::string("action ") + action_to_string(a) + " " + action_label(a));

        // arguments
        desc += "(";
        bool need_comma = false;
        for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
            if( model_[args(a, mo)] ) {
                if( need_comma ) desc += ",";
                desc += meta_object_to_string(mo);
                need_comma = true;
            }
        }
        desc += "):";

        // static precondition
        desc += " static-pre={";
        need_comma = false;
        for( Unary u = Unary(0); u < num_static_unary_predicates_; ++u ) {
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
                if( model_[unary(u, a, mo)] ) {
                    if( need_comma ) desc += ",";
                    desc += unary_schema_to_string(u, mo);
                    need_comma = true;
                }
            }
        }
        for( Binary b = Binary(0); b < num_static_binary_predicates_; ++b ) {
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
                for( MetaObject mop = MetaObject(0); mop < num_meta_objects_; ++mop ) {
                    if( model_[binary(b, a, mo, mop)] ) {
                        if( need_comma ) desc += ",";
                        desc += binary_schema_to_string(b, mo, mop);
                        need_comma = true;
                    }
                }
            }
        }
        desc += "}";

        // precondition
        desc += ", pre={";
        need_comma = false;
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            if( model_[pre0(a, mk)] || model_[pre1(a, mk)] ) {
                if( need_comma ) desc += ",";
                if( model_[pre0(a, mk)] ) desc += "-";
                desc += atom_schema_desc(mk);
                need_comma = true;
            }
        }
        desc += "}";

        // effect
        desc += ", eff={";
        need_comma = false;
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            if( model_[eff0(a, mk)] || model_[eff1(a, mk)] ) {
                if( need_comma ) desc += ",";
                if( model_[eff0(a, mk)] ) desc += "-";
                desc += atom_schema_desc(mk);
                need_comma = true;
            }
        }
        desc += "}";

        return desc;
    }

    // description of ground action
    std::string ground_action_desc(const std::vector<int> &tuple) const {
        assert(1 + num_meta_objects_ == int(tuple.size()));
        Action a = Action(tuple.front());
        std::string desc(action_label(a) + "(");
        bool need_comma = false;
        for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
            if( model_.at(args(a, mo)) ) {
                if( need_comma ) desc += ",";
                desc += object_to_string(Object(tuple.at(1 + mo)));
                need_comma = true;
            }
        }
        desc += ")";
        return desc;
    }
    std::string ground_action_desc(Transition t) const {
        std::string desc;
        bool found = false;
        auto foo = [this, t, &desc, &found](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            if( model_.at(G(tuple)) && (tuple.at(0) == t) ) {
                assert(!found);
                std::vector<int> rtuple(&tuple[1], &tuple[2 + num_meta_objects_]);
                desc += ground_action_desc(rtuple);
                found = true;
            }
        };
        G.enumerate_vars_from_multipliers(foo);
        return desc;
    }
    void ground_action_desc(std::ostream &os, Layer l, Action a) const {
        auto foo = [this, &os, l, a](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
            if( (tuple.at(0) == l) && (tuple.at(1) == a) ) {
                if( model_.at(gtuple(tuple)) ) {
                    std::vector<int> ntuple { a };
                    ntuple.insert(ntuple.end(), &tuple[2], &tuple[tuple.size()]);
                    os << ground_action_desc(ntuple) << std::endl;
                }
            }
        };
        gtuple.enumerate_vars_from_multipliers(foo);
    }

    // state to string
    std::string state_to_string(State s) const {
        return std::string("s") + std::to_string(s);
    }

    // transition to string
    std::string transition_to_string(Transition t) const {
        return std::string("t") + std::to_string(t);
    }

    // description of state
    std::string state_desc(State s, bool atoms_desc = true, bool appl_desc = false) const {
        Layer l = s_layer_.at(s);
        std::string desc(state_to_string(s));
        if( atoms_desc ) {
            desc += "{";
            bool need_comma = false;
            for( int i = 0; i < int(features_per_layer_[l].size()); ++i ) {
                Feature k = features_per_layer_[l].at(i);
                if( non_constant_feature(k) && model_.at(phi(k, s)) ) {
                    if( need_comma ) desc += ",";
                    desc += ground_atom_desc(k);
                    need_comma = true;
                }
            }
            desc += "}";
        } else {
            desc += "/";
            for( int i = 0; i < int(features_per_layer_[l].size()); ++i ) {
                Feature k = features_per_layer_[l].at(i);
                if( non_constant_feature(k) )
                    desc += std::to_string(model_.at(phi(k, s)));
            }
        }

        // applicability of actions
        if( appl_desc ) {
            desc += ": appl={";
            if( encoding("applicable-actions-alt") ) {
                bool need_comma = false;
                for( Action a = Action(0); a < num_actions_; ++a ) {
                    for( int i = 0; i < int(transitions_per_layer_[l].size()); ++i ) {
                        Transition t = transitions_per_layer_[l].at(i);
                        if( model_.at(appl(a, t, s)) ) {
                            if( need_comma ) desc += ",";
                            desc += action_to_string(a) + "/" + transition_to_string(t);
                            need_comma = true;
                        }
                    }
                }
            } else {
                assert(encoding("applicable-actions-tuples"));

                // output a(<tuple>)/<transition>
                bool need_comma = false;
                auto foo = [this, &desc, &need_comma, s](const SAT::VarSet &varset, const std::vector<int> &tuple) -> void {
                    State sp(tuple.at(1 + num_meta_objects_));
                    std::vector<int> rtuple(&tuple[0], &tuple[1 + num_meta_objects_]);
                    if( model_.at(appl2(tuple)) && (sp == s) ) {
                        if( need_comma ) desc += ",";
                        desc += ground_action_desc(rtuple);
                        desc += "/";
                        bool need_comma2 = false;
                        auto bar = [this, &desc, &need_comma2, s, &rtuple](const SAT::VarSet &varset, const std::vector<int> &otuple) -> void {
                            Transition t(otuple.at(0));
                            std::vector<int> rotuple(&otuple[1], &otuple[2 + num_meta_objects_]);
                            if( model_.at(G(otuple)) && (tr_src_.at(t) == s) && (rtuple == rotuple) ) {
                                if( need_comma2 ) desc += ",";
                                desc += transition_to_string(t); // + "/" + ground_action_desc(rotuple);
                                need_comma2 = true;
                            }
                        };
                        G.enumerate_vars_from_multipliers(bar); // G(t,a,<tuple>)
                        need_comma = true;
                    }
                };
                appl2.enumerate_vars_from_multipliers(foo); // appl(a,<tuple>,s)
            }
            desc += "}";
        }

        return desc;
    }


    // description of transitions
    std::string transition_desc(Transition t, bool mapt_desc = false, bool mapf_desc = false) const {
        Layer l = tr_layer_.at(t);
        State src = tr_src_.at(t);
        State dst = tr_dst_.at(t);
        std::string desc(transition_to_string(t) + ":");
        desc += " src=" + state_desc(src);
        desc += ", dst=" + state_desc(dst);

        // ground action
        desc += std::string(", ground-action=") + ground_action_desc(t);
        Action action = Action(-1);
        for( Action a = Action(0); a < num_actions_; ++a ) {
            if( model_.at(map(t, a)) ) {
                action = a;
                break;
            }
        }
        assert(action != -1);
        desc += ", action=" + action_to_string(action);

        // mapf description
        if( mapf_desc ) {
            desc += ", mapf={";
            bool need_comma = false;
            for( int i = 0; i < int(features_per_layer_[l].size()); ++i ) {
                Feature k = features_per_layer_[l].at(i);
                for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
                    if( model_.at(mapf(t, k, mk)) ) {
                        if( need_comma ) desc += ",";
                        desc += feature_to_string(k) + "->" + meta_feature_to_string(mk);
                        need_comma = true;
                    }
                }
            }
            desc += "}";
        }

        // mapt description
        if( mapt_desc ) {
            desc += ", mapt={";
            bool need_comma = false;
            for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo ) {
                for( Object o = Object(0); o < objects_per_layer_[l]; ++o ) {
                    if( model_.at(mapt(t, mo, o)) ) {
                        if( need_comma ) desc += ",";
                        desc += meta_object_to_string(mo) + "->" + object_to_string(o);
                        need_comma = true;
                    }
                }
            }
            desc += "}";
        }

        return desc;
    }

    void decode_meta_layer(std::ostream &os) const {
        assert(satisfiable_ && (model_.size() == variables_.size()));

        // output instantiation for atoms that define meta-layer
        std::vector<const SAT::VarSet*> meta_layer {
            &arity, &using1, &using2, &label,
            &pre0, &pre1, &eff0, &eff1,
            &args, &atom1, &atom2,
            &unary, &binary
        };
        for( size_t i = 0; i < meta_layer.size(); ++i )
            meta_layer[i]->print(os, model_, true);
    }

    void decode_model(std::ostream &os) const override {
        assert(satisfiable_ && (model_.size() == variables_.size()));

        os << std::endl << red("Meta-layer:", !options_.disable_colors_) << std::endl;

        // meta-features
        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk )
            os << meta_feature_desc(mk) << std::endl;

        // actions in terms of meta-features
        os << std::endl << red("Action in terms of meta-features:", !options_.disable_colors_) << std::endl;
        for( Action a = Action(0); a < num_actions_; ++a )
            os << action_mk_desc(a) << std::endl;

        // action schema
        os << std::endl << red("Action schema:", !options_.disable_colors_) << std::endl;
        for( Action a = Action(0); a < num_actions_; ++a )
            os << action_schema_desc(a) << std::endl;

        // static predicates
        os << std::endl << red("Extension of unary static predicates:", !options_.disable_colors_) << std::endl;
        r.print(os, model_);

        os << std::endl << red("Extension of binary static predicates:", !options_.disable_colors_) << std::endl;
        s.print(os, model_);

#if 0
        // ground actions
        os << std::endl << red("Ground actions:", !options_.disable_colors_) << std::endl;
        gtuple.print(os, model_);

        // meta-feature vectors (lex ordering)
        os << std::endl;
        os << "lexical ordering of meta-features: using(mk)";
        if( options_.disjoint_meta_features_ ) os << " using(a,mk)";
        os << " atom(mk,p) atom(mk,i,mo)";
        if( !options_.disjoint_meta_features_ ) os << " using(a,mk)";
        os << std::endl;

        for( MetaFeature mk = MetaFeature(0); mk < num_meta_features_; ++mk ) {
            os << "mk" << mk << ": ";

            // used by some action
            os << model_.at(using1(mk)) << " ";

            if( options_.disjoint_meta_features_ ) {
                // usage in actions
                for( Action a = Action(0); a < num_actions_; ++a )
                    os << model_.at(using2(a, mk));
                os << " ";
            }

            // assignment to atoms
            for( Atom p = Atom(0); p < num_atoms_; ++p )
                os << model_.at(atom1(mk, p));
            os << " ";

            // assignment of meta-objects to arguments
            for( Arity i = Arity(0); i < 1 + max_arity_; ++i ) {
                for( MetaObject mo = MetaObject(0); mo < num_meta_objects_; ++mo )
                    os << model_.at(atom2(mk, i, mo));
                os << " ";
            }

            if( !options_.disjoint_meta_features_ ) {
                // usage in actions
                for( Action a = Action(0); a < num_actions_; ++a )
                    os << model_.at(using2(a, mk));
            }

            os << std::endl;
        }
        os << std::endl;
#endif

        // layers
        for( Layer l = Layer(0); l < num_layers_; ++l ) {
            os << std::endl << red(std::string("Layer ") + std::to_string(l) + ":", !options_.disable_colors_) << std::endl;

            // features
            for( int i = 0; i < int(features_per_layer_[l].size()); ++i ) {
                Feature k = Feature(features_per_layer_[l].at(i));
                if( non_constant_feature(k) )
                    os << feature_desc(k) << std::endl;
            }

            // ground actions
            os << std::endl << red("Ground actions:", !options_.disable_colors_) << std::endl;
            for( Action a = Action(0); a < num_actions_; ++a )
                ground_action_desc(os, l, a);

            // States + Appl
            os << std::endl << red("States:", !options_.disable_colors_) << std::endl;
            for( int i = 0; i < int(states_per_layer_[l].size()); ++i )
                os << state_desc(states_per_layer_[l].at(i), true, true) << std::endl;

            // transitions
            os << std::endl << red("Transitions:", !options_.disable_colors_) << std::endl;
            for( int i = 0; i < int(transitions_per_layer_[l].size()); ++i )
                os << transition_desc(transitions_per_layer_[l].at(i), false, false) << std::endl;

#if 0
            // calculate set of affected features per action label
            std::vector<std::set<std::set<Feature> > > affected_features(num_actions_);
            std::vector<Label> action_labels(num_actions_);
            for( Action a = Action(0); a < num_actions_; ++a ) {
                for( Label l = Label(0); l < num_labels_; ++l ) {
                    if( model_.at(label(a, l)) )
                        action_labels[a] = l;
                }

                for( int i = 0; i < int(transitions_per_layer_[l].size()); ++i ) {
                    Transition t = transitions_per_layer_[l].at(i);
                    State src = tr_src_[t];
                    State dst = tr_dst_[t];
                    if( model_.at(map(t, a)) ) {
                        std::set<Feature> features;
                        for( int j = 0; j < int(features_per_layer_[l].size()); ++j ) {
                            Feature k = features_per_layer_[l].at(j);
                            if( model_.at(phi(k, src)) != model_.at(phi(k, dst)) ) {
                                features.insert(k);
                            }
                        }
                        affected_features[a].insert(features);
                    }
                }
            }

            // instantiation of action schemas in layer
            for( Action a = Action(0); a < num_actions_; ++a ) {
                os << "label " << label_as_string_[action_labels[a]]
                   << ": action=a" << a
                   << ", aff-features={";

                bool need_comma = false;
                for( std::set<std::set<Feature> >::const_iterator it = affected_features[a].begin(); it != affected_features[a].end(); ++it ) {
                    if( need_comma ) os << ",";
                    os << " {";
                    bool need_comma2 = false;
                    for( std::set<Feature>::const_iterator jt = it->begin(); jt != it->end(); ++jt ) {
                        if( need_comma2 ) os << ",";
                        os << feature_to_string(*jt);
                        need_comma2 = true;
                    }
                    os << "}";
                    need_comma = true;
                }
                need_comma = false;

                os << " }, trans={";
                for( int i = 0; i < int(transitions_per_layer_[l].size()); ++i ) {
                    Transition t = transitions_per_layer_[l].at(i);
                    if( model_.at(map(t, a)) ) {
                        if( need_comma ) os << ",";
                        os << transition_to_string(t);
                        need_comma = true;
                    }
                }
                need_comma = false;
                os << "}" << std::endl;
            }
            os << std::endl;
#endif
        }
    }
};

}; // namespace

#endif

