#ifndef POSETOFSETS_HPP_INCLUDED
#define POSETOFSETS_HPP_INCLUDED

#include <iostream>
#include <boost/lexical_cast.hpp>
#define prt(x) std::cerr << #x " = '" << x << "'" << std::endl;
#define prt2(x,y) std::cerr << #x " = '" << x << "'\t" << #y " = '" << y << "'" << std::endl;
#define prt3(x,y,z) std::cerr << #x " = '" << x << "'\t" << #y " = '" << y << "'\t" << #z " = '" << z << "'" << std::endl;
#define prt4(x,y,z,u) std::cerr << #x " = '" << x << "'\t" << #y " = '" << y << "'\t" << #z " = '" << z << "'\t" << #u " = '" << u << "'" << std::endl;
#define prt5(x,y,z,u,v) std::cerr << #x " = '" << x \
    << "'\t" << #y " = '" << y \
    << "'\t" << #z " = '" << z \
    << "'\t" << #u " = '" << u \
    << "'\t" << #v " = '" << v \
    << "'" << std::endl;
#define prt6(x,y,z,u,v,w) std::cerr << #x " = '" << x \
    << "'\t" << #y " = '" << y \
    << "'\t" << #z " = '" << z \
    << "'\t" << #u " = '" << u \
    << "'\t" << #v " = '" << v \
    << "'\t" << #w " = '" << w \
    << "'" << std::endl;
#define prt7(x,y,z,u,v,w,t) std::cerr << #x " = '" << x \
    << "'\t" << #y " = '" << y \
    << "'\t" << #z " = '" << z \
    << "'\t" << #u " = '" << u \
    << "'\t" << #v " = '" << v \
    << "'\t" << #w " = '" << w \
    << "'\t" << #t " = '" << t \
    << "'" << std::endl;
#define pt(x)  "\"" << #x << "\" = " << x << ";\t"

#define assert2(expr, value)\
   if (!(expr)){\
        if ("\"" + boost::lexical_cast<std::string>(value) + "\"" == __STRING(value) ){\
            std::cerr << (value) << std::endl;\
        }else{\
            std::cerr << (__STRING(value)) << ": " << (value) << std::endl;\
        }\
        __assert_fail (__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION); abort(); \
   }

#define assertss(expr, value)\
   if (!(expr)){\
       std::cerr << "BOEM!" << std::endl << value << std::endl; \
        __assert_fail (__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION); abort(); \
   }
#define assertssexec(expr, value, __e__)\
   if (!(expr)){\
       std::cerr << "BOEM!" << std::endl << value << std::endl; __e__; \
        __assert_fail (__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION); abort(); \
   }

#define assert3(expr, v1, v2)\
   if (!(expr)){\
        if ("\"" + boost::lexical_cast<std::string>(v1) + "\"" == __STRING(v1) ){\
            std::cerr << (v1) << std::endl;\
        }else{\
            std::cerr << (__STRING(v1)) << ": " << (v1) << std::endl;\
        }\
        if ("\"" + boost::lexical_cast<std::string>(v2) + "\"" == __STRING(v1) ){\
            std::cerr << (v2) << std::endl;\
        }else{\
            std::cerr << (__STRING(v2)) << ": " << (v2) << std::endl;\
        }\
        __assert_fail (__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION);  abort();\
   }

template<typename _InputIterator1, typename _InputIterator2>
bool has_intersection(_InputIterator1 __first1, _InputIterator1 __last1, _InputIterator2 __first2, _InputIterator2 __last2)
{
    while(__first1 != __last1 && __first2 != __last2)
        if(*__first1 < *__first2)
            ++__first1;
        else if(*__first2 < *__first1)
            ++__first2;
        else
            return true;

    return false;
}

#include <vector>
#include <map>
#include <limits>
#include <string>
#include <assert.h>
#include <iostream>
#include <set>
#include <list>
#include <functional>
#include <sstream>
#include <algorithm>
#include <stack>

#include "prettyprint.hpp"

void writeToFile(const std::string file, const std::string& s);

namespace
{
    using std::vector;
    using std::size_t;
    using std::map;
    using std::set;
    using std::list;

    template<typename traits>
    class PosetOfSets;
    template<typename traits>
    class PosetNode;
    template<typename traits>
    class PosetQuery;

    template<typename traits>
    class PosetOfSets
    {
    public:
        PosetNode<traits> all, //represents the ideal node with size infinity
                empty; //the opposite, size 0.

        std::set<PosetNode<traits>*> candidate_nodes_for_deletion;
        std::multimap<typename traits::subelement_type, PosetNode<traits>*> reverse_introductions;

        typedef PosetNode<traits> node;

    public:
        PosetOfSets();
        ~PosetOfSets();

        PosetNode<traits>*
        getNodeCreateIfMissing2(typename traits::subelementlookup& subelements, const bool introduce_common_subsets = true);

        PosetNode<traits>* getAll() const;
        PosetNode<traits>* getEmpty() const;

        typename traits::value_type* getValue(PosetNode<traits>* node);
        void setValue(PosetNode<traits>* node, typename traits::value_type value);

        void discardNode(PosetNode<traits>* node);

        std::string toDot() const;

        void dump();

    };

    template<typename traits>
    void PosetOfSets<traits>::dump()
    {
        writeToFile("poset.dot", toDot());
    }

    template<typename traits>
    class PosetNode
    {
    public:
        typedef PosetNode<traits> node;

        // these define the DAG:
        std::vector<node*> higher, lower;

        // these define the spanning tree:
        std::vector<node*> deriving_chains;
        node* lower_chain;

        std::set<typename traits::subelement_type> elements_introduced_here_relative_to_lower_chain;
        const size_t size;

        // the payload:
        typename traits::value_type value;

        PosetNode(const size_t size_, const PosetOfSets<traits> * myposs, const bool attached_);

        void applyDownwardsInChain(PosetOfSets<traits> * myposs, typename std::function<bool(PosetNode<traits>*)> f);

        void cascadingDelete(size_t chain);

    public:

        size_t getSize() const;

        std::string toDot() const
        {
            std::stringstream oss;
            oss //
            << "this: " << this << "\\\l" //
                    << "larger: " << higher << "\\\l" //
                    << "deriving_chains: " << deriving_chains << "\\\l" //
                    << "size: " << size << "\\\l" //
                    << "smaller: " << "lower_chain: " << lower_chain << "\\\l" //
                    << "introduced here: " << elements_introduced_here_relative_to_lower_chain << "\\\l" //
                    << lower << "\\\l";

//                    << pt(size) << pt(my_chain) << pt(lower_chain) << pt(higher_chain) << pt(higher) << pt(converging_chains) << pt(deriving_chains) << pt(lower);
            return oss.str();
        }

        std::string toString() const
        {
            std::stringstream oss;
            oss << pt(higher) << std::endl //
                    << pt(deriving_chains) << std::endl //
                    << pt(size) << std::endl //
                    << pt(elements_introduced_here_relative_to_lower_chain) << std::endl //
                    << pt(lower_chain) << std::endl //
                    << pt(lower) << std::endl;
            return oss.str();

        }
    };

    template<typename traits>
    class PosetQuery
    {
    public:

//        friend class PosetOfSets<traits> ;
//        friend class vector<PosetQuery<traits> > ;

        PosetOfSets<traits> * myposs;

        typedef PosetNode<traits> node;

//        std::map<std::size_t, node*> upperbounds_to_be_lowered, lowerbounds_to_be_lifted_ordered;
    public:

        PosetQuery(PosetOfSets<traits> * mp) :
                myposs(mp)
        {
        }

        void find_smallest_sets_surrounding_given_set(const std::set<typename traits::subelement_type>& subelements, std::set<node*>& current_intersection)
        {
            // it looks ok but ... should i refer to the somewhat larger set when it is from a different chain too ?
            if(subelements.empty())
            {
                current_intersection.insert(myposs->getEmpty());
                return;
            }

            for(auto& se : subelements)
            {
                auto iterpair = myposs->reverse_introductions.equal_range(se);
                if(iterpair.first == iterpair.second)
                {
                    // std::cerr << "wth, no subelement " << se << std::endl;
                    current_intersection.clear();
                    current_intersection.insert(myposs->getAll());
                    return;
                }
                if(current_intersection.empty())
                { // initial element ...
                    for(auto j = iterpair.first; j != iterpair.second; ++j)
                    {
                        current_intersection.insert(j->second);
                    }
                }
                else
                {
                    std::set<node*> next_intersection;
                    for(auto j = iterpair.first; j != iterpair.second; ++j)
                    {
                        for(node* k : current_intersection)
                        {
                            node* largenode = k;
                            node * smallnode = j->second;
                            if(k->getSize() < smallnode->getSize())
                            {
                                std::swap(largenode, smallnode);
                            }
                            bool found_the_small_node_beneath_the_large_node = false;
                            largenode->applyDownwardsInChain(myposs, [&](node* ancestor)
                            {
                                if (ancestor == smallnode)
                                {
                                    found_the_small_node_beneath_the_large_node = true;
                                }
                                return !found_the_small_node_beneath_the_large_node;
                            });
                            if(found_the_small_node_beneath_the_large_node)
                            {
                                next_intersection.insert(largenode);
                                // do not break here, there could be more ... larger nodes ...
                            }
                        }
                    }
                    current_intersection.clear();
                    if(next_intersection.empty())
                    {
                        current_intersection.insert(myposs->getAll());
                        return;
                    }
                    else
                    {
                        current_intersection.insert(next_intersection.begin(), next_intersection.end());
                    }
                }
            }

            // THESE ARE THE PER QUEUE SMALLEST SURROUNDING!
            // still got to get rid of a few!

            std::map<std::size_t, std::set<node*> > size_to_candidate_nodes;
            for(auto i : current_intersection)
            {
//                again.insert(i);
                size_to_candidate_nodes[i->getSize()].insert(i);
            }
            for(auto i = size_to_candidate_nodes.begin(); i != size_to_candidate_nodes.end(); ++i)
            {
                for(node* rootnode : i->second)
                {
                    std::stack<node*> todo;
                    todo.push(rootnode);

                    while(!todo.empty())
                    {
                        node* current_ok_node = todo.top();
                        todo.pop();
                        for(node* k : current_ok_node->higher)
                        {
                            //prt(k.second->getSize());
                            assert(k->getSize() > rootnode->getSize());
                            auto zucht = size_to_candidate_nodes.find(k->getSize());
                            if(zucht != size_to_candidate_nodes.end())
                            {
                                //prt2(zucht->second, k.second);
                                zucht->second.erase(k);
                                if(zucht->second.empty())
                                {
                                    //std::cerr << "higher set for size " << zucht->first << " is now empty" << std::endl;
                                    size_to_candidate_nodes.erase(zucht);
                                }
                            }
                            if(k->getSize() <= size_to_candidate_nodes.rbegin()->first)
                            {
                                todo.push(k);
                            }
                        }
                    }
                }
            }

            current_intersection.clear();
            for(auto &i : size_to_candidate_nodes)
            {
                current_intersection.insert(i.second.begin(), i.second.end());
            }

        }

        void find_largest_subsets_of_given_set(const std::set<typename traits::subelement_type>& subelements, std::set<node*>& ret)
        {
            if(!find_largest_subsets_of_given_set_(subelements, myposs->getEmpty(), ret))
            {
                //  std::cerr << "no largestsubsets??" << std::endl;
                ret.insert(myposs->getEmpty());
            }

            // these are the largest subsets PER QUEUE ... there is still some filtering to be done ...

            std::map<std::size_t, std::set<node*> > size_to_candidate_nodes;
            for(node* i : ret)
            {
                size_to_candidate_nodes[i->getSize()].insert(i);
            }
            for(auto i = size_to_candidate_nodes.rbegin(); i != size_to_candidate_nodes.rend(); ++i)
            {
                for(node* rootnode : i->second)
                {
                    std::stack<node*> todo;
                    todo.push(rootnode);

                    while(!todo.empty())
                    {
                        node* current_ok_node = todo.top();
                        todo.pop();
                        for(node* k : current_ok_node->lower)
                        {
                            //prt(k.second->getSize());
                            assert(k->getSize() < rootnode->getSize());
                            auto zucht = size_to_candidate_nodes.find(k->getSize());
                            if(zucht != size_to_candidate_nodes.end())
                            {
                                //prt2(zucht->second, k.second);
                                zucht->second.erase(k);
                                if(zucht->second.empty())
                                {
                                    //std::cerr << "higher set for size " << zucht->first << " is now empty" << std::endl;
                                    size_to_candidate_nodes.erase(zucht);
                                }
                            }
                            if(k->getSize() >= size_to_candidate_nodes.begin()->first)
                            {
                                todo.push(k);
                            }
                        }
                    }
                }
            }

            ret.clear();
            for(auto &i : size_to_candidate_nodes)
            {
                ret.insert(i.second.begin(), i.second.end());
            }

        }

        bool find_largest_subsets_of_given_set_(const std::set<typename traits::subelement_type>& subelements, node* current, std::set<node*>& ret)
        {
            // there wont be a result for every chain
            const auto& entries = current->elements_introduced_here_relative_to_lower_chain;
            bool i_can_replace_below = entries.empty() || std::includes(subelements.begin(), subelements.end(), entries.begin(), entries.end());
            if(i_can_replace_below)
            {
                bool i_could_be_replaced = false;
                for(auto& higher_node : current->deriving_chains)
                {
                    if(higher_node->getSize() <= subelements.size())
                    {
                        i_could_be_replaced |= find_largest_subsets_of_given_set_(subelements, higher_node, ret);
                    }
                }
//                prt4(i_can_replace_below, i_could_be_replaced, current, myposs->getEmpty());
                if(!i_could_be_replaced)
                {
                    //std::cerr << "Largest subset (for this part of the spanning tree) is " << current->toString() << pt(subelements) << std::endl;
                    ret.insert(current);
                }
            }
            return i_can_replace_below;
        }

    public:
        PosetQuery(const PosetQuery<traits>& pca);
    protected:

        bool isComparable(const size_t chain) const;
        bool couldChangeStill(const size_t chain) const;
        PosetNode<traits>* couldItBeTheLowerBound(const size_t size) const;
        void getLowerAndUpper(PosetNode<traits> * & l, PosetNode<traits> * & u) const;
    };

    template<typename traits>
    PosetOfSets<traits>::PosetOfSets() :
            all(std::numeric_limits<size_t>::max(), this, true),
            empty(0, this, true)
    {
        all.lower.push_back(&empty);
        empty.higher.push_back(&all);
        // do not add "All" to the the spanning tree
    }

    template<typename traits>
    PosetOfSets<traits>::~PosetOfSets()
    {
        // what about using unique_ptr's for the spanning tree?
        std::set<node*> all_nodes;
        for(auto& i : reverse_introductions)
        {
            all_nodes.insert(i.second);
        }
        for(node* n : all_nodes)
        {
            delete n;
        }
    }

    template<typename _InputIterator1, typename _InputIterator2>
    void iterate_intersection(_InputIterator1 __first1, _InputIterator1 __last1, _InputIterator2 __first2, _InputIterator2 __last2, std::function<bool(_InputIterator1, _InputIterator2)>& f)
    {
        while(__first1 != __last1 && __first2 != __last2)
            if(__first1->first < __first2->first)
            {
                ++__first1;
            }
            else
            {
                if(!(__first2->first < __first1->first))
                {
                    if(!f(__first1->second, __first2->second))
                    {
                        return;
                    }
                }
                ++__first2;
            }

        return false;
    }

    template<typename T>
    void eraseValueFromVector(std::vector<T>& c, T t)
    {
        for(auto i = c.begin(); i != c.end();)
        {
            if(*i == t)
            {
                i = c.erase(i);
            }
            else
            {
                ++i;
            }
        }
    }

    template<typename T>
    bool containedInVector(std::vector<T>& c, T t)
    {
        for(auto i = c.begin(); i != c.end(); ++i)
        {
            if(*i == t)
            {
                return true;
            }
        }
        return false;
    }

    template<typename traits>
    PosetNode<traits>*
    PosetOfSets<traits>::getNodeCreateIfMissing2(typename traits::subelementlookup& subelements, const bool introduce_common_subsets)
    {
//        prt(subelements);
        PosetQuery<traits> query(this);
        std::set<node*> largest_subsets, smallest_sets_surrounding;
        query.find_largest_subsets_of_given_set(subelements, largest_subsets);
        for(auto& subset : largest_subsets)
        {
            if(subset->getSize() == subelements.size())
            {
                return subset;
            }
            assert(subelements.size() >= subset->getSize());
        }

        query.find_smallest_sets_surrounding_given_set(subelements, smallest_sets_surrounding);
        for(auto& surrounding : smallest_sets_surrounding)
        {
            if(surrounding->getSize() == subelements.size())
            {
                //existing_solution_node = surrounding.second;
                prt3(subelements, largest_subsets, smallest_sets_surrounding);
                assertssexec(false, pt(subelements) << std::endl << pt(surrounding) << std::endl << pt(smallest_sets_surrounding) << std::endl << pt(largest_subsets), dump()); // can't happen ... should have fired in the previous test ...
                return surrounding;
            }
            assertssexec(surrounding->getSize() > subelements.size(), pt(surrounding->toDot()) << std::endl << pt(subelements), dump());
        }

        assert(!largest_subsets.empty());
        assert(!smallest_sets_surrounding.empty());

        node* n = new node(subelements.size(), this, true);

        node* beste_parent = 0;

        for(node* candidate_chain_bringer : largest_subsets) // kies de kandidaatcluster zo dat we zo weinig mogelijk subelementen moeten introduceren
        {
            if(beste_parent == 0 || candidate_chain_bringer->getSize() > beste_parent->getSize())
            {
                beste_parent = candidate_chain_bringer;
            }
            else if(candidate_chain_bringer->getSize() == beste_parent->getSize() && candidate_chain_bringer->deriving_chains.size() > beste_parent->deriving_chains.size())
            {
                // the benefit of taking from a parent with derived chains ... is that it is cheaper to build the introduced set...
                beste_parent = candidate_chain_bringer;
            }
        }
        assert(beste_parent); // i might get into his chain or i might derive from him...
        n->lower_chain = beste_parent;

        for(node* node_deriving_from_beste_parent : beste_parent->deriving_chains)
        {
//            node* node_deriving_from_beste_parent = beste_parent->getHigherFromNode(chain_deriving_from_beste_parent);
            if(node_deriving_from_beste_parent != getAll())
            {
                if(smallest_sets_surrounding.find(node_deriving_from_beste_parent) != smallest_sets_surrounding.end())
                {
                    // use this node, which might be completely unrelated to the rest of the story, to build up our introduced set
                    for(auto introduced : node_deriving_from_beste_parent->elements_introduced_here_relative_to_lower_chain)
                    {
                        if(subelements.find(introduced) != subelements.end())
                        {
                            n->elements_introduced_here_relative_to_lower_chain.insert(introduced);
                        }
                    }
                    assertssexec(!n->elements_introduced_here_relative_to_lower_chain.empty(), pt(subelements), dump());
                    break;
                }
            }
        }

        if(n->elements_introduced_here_relative_to_lower_chain.empty())
        {
            /// we could not find a donor, no node could tell us what should be in the introduced set ... so do it the hard way ...

            n->elements_introduced_here_relative_to_lower_chain.insert(subelements.begin(), subelements.end());

            //std::size_t total = 0;
            beste_parent->applyDownwardsInChain(this, [&](PosetNode<traits>* ancestor)
            {
                //prt2(ancestor,ancestor->elements_introduced_here_relative_to_lower_chain);
                //total += ancestor->elements_introduced_here_relative_to_lower_chain.size();
                    assertss(std::includes(n->elements_introduced_here_relative_to_lower_chain.begin(),n->elements_introduced_here_relative_to_lower_chain.end(), ancestor->elements_introduced_here_relative_to_lower_chain.begin(), ancestor->elements_introduced_here_relative_to_lower_chain.end()),pt(n->elements_introduced_here_relative_to_lower_chain) << std::endl << pt(ancestor->elements_introduced_here_relative_to_lower_chain));
                    for (auto& j : ancestor->elements_introduced_here_relative_to_lower_chain)
                    {
                        n->elements_introduced_here_relative_to_lower_chain.erase(j);
                    }
                    //                                n->elements_introduced_here_relative_to_lower_chain.erase(ancestor->elements_introduced_here_relative_to_lower_chain.begin(), ancestor->elements_introduced_here_relative_to_lower_chain.end());

                    return true;// keep going
                });

        }
        assert(!n->elements_introduced_here_relative_to_lower_chain.empty());

        for(node* surrounding : smallest_sets_surrounding)
        {
            if(surrounding != getAll())
            {
                n->deriving_chains.push_back(surrounding);

                for(auto i = surrounding->elements_introduced_here_relative_to_lower_chain.begin(); i != surrounding->elements_introduced_here_relative_to_lower_chain.end();)
                {
                    if(subelements.find(*i) != subelements.end())
                    {
                        auto iterpair = reverse_introductions.equal_range(*i);
                        //bool removed_reverse_link = false;
                        for(auto& j = iterpair.first; j != iterpair.second; ++j)
                        {
                            if(j->second == surrounding)
                            {
                                j = reverse_introductions.erase(j);
                                break;
                            }
                        }
                        i = surrounding->elements_introduced_here_relative_to_lower_chain.erase(i);
                    }
                    else
                    {
                        ++i;
                    }
                }

                while(surrounding->getSize() > surrounding->elements_introduced_here_relative_to_lower_chain.size() + n->getSize())
                {
                    // go find the missing elements in the basement of surrounding...// TODO: WE REALLY REALLY NEED A NEW CHAININDEX IF THIS HAPPENS

                    surrounding->lower_chain->applyDownwardsInChain(this, [&](node* ancestor)
                    {
                        for (auto& grr : ancestor->elements_introduced_here_relative_to_lower_chain)
                        {
                            if (subelements.find(grr) == subelements.end())
                            {
                                //          std::cerr << "what makes me special: "; prt(grr);
                            surrounding->elements_introduced_here_relative_to_lower_chain.insert(grr);
                            //todo: modify the reverse connections but .. with the new chainindex (which i dont have yet)
                            reverse_introductions.insert(std::make_pair(grr, surrounding ));
                        }
                    }

                    return surrounding->getSize() > surrounding->elements_introduced_here_relative_to_lower_chain.size() + n->getSize();});
                }

                assertssexec(!surrounding->elements_introduced_here_relative_to_lower_chain.empty(), pt(surrounding->getSize()) << std::endl << pt(subelements)<< std::endl << pt(surrounding->lower),
                {
                    reverse_introductions.insert(std::make_pair(0, surrounding))
                    ;
                    dump()
                    ;
                });

                eraseValueFromVector(surrounding->lower_chain->deriving_chains, surrounding);
                surrounding->lower_chain = n;
            }
        }

        // remove all ... direct connections ... between the higher and lower sets ... because the new node is going to take over ...
        // ---> cutting chains here ... how to fix them???
        for(node* above : smallest_sets_surrounding)
        {
            for(node* below : largest_subsets)
            {
                eraseValueFromVector(above->lower, below);
                eraseValueFromVector(below->higher, above);
            }
        }

        // THE REAL CONNECTING:
//        prt3(largest_subsets, smallest_sets_surrounding, n);
        for(node* below : largest_subsets)
        {
            below->higher.push_back(n);
            n->lower.push_back(below);
        }

        for(node* above : smallest_sets_surrounding)
        {
            above->lower.push_back(n);
            n->higher.push_back(above);
        }

        for(auto& i : n->elements_introduced_here_relative_to_lower_chain)
        {
            reverse_introductions.insert(std::make_pair(i, n));
        }
        beste_parent->deriving_chains.push_back(n);

        // this new node ... the parent of it ... ...
        // does it have multiple derived chains?
        // can i reduce the amount of elements in their introduced sets?
        if(introduce_common_subsets && largest_subsets.size() == 1)
        {
            std::set<node*> siblings;
            for(auto i : beste_parent->deriving_chains)
            {
                siblings.insert(i);
            }
            if(siblings.size() > 1)
            {
                node* a = 0;
                node* b = 0;
                for(auto i : siblings)
                {
                    if(a == 0 || i->getSize() >= a->getSize())
                    {
                        a = i;
                    }
                }
                for(auto i : siblings)
                {
                    if(i != a && (b == 0 || i->getSize() >= b->getSize()))
                    {
                        b = i;
                    }
                }
                assertss(a != b, pt(siblings) << std::endl //
                        << pt(a) << std::endl//
                        << pt(b) << std::endl//
                        << pt(a->getSize()) << std::endl//
                        << pt(b->getSize()));
                assert(beste_parent);

                std::size_t numa = a->elements_introduced_here_relative_to_lower_chain.size();
                std::size_t numb = b->elements_introduced_here_relative_to_lower_chain.size();

                std::set<typename traits::subelement_type> intersection;
                std::set_intersection(a->elements_introduced_here_relative_to_lower_chain.begin(), a->elements_introduced_here_relative_to_lower_chain.end(), //
                        b->elements_introduced_here_relative_to_lower_chain.begin(), b->elements_introduced_here_relative_to_lower_chain.end(), //
                        std::inserter(intersection, intersection.begin()));

                if(!intersection.empty() && intersection.size() < numa && intersection.size() < numb) //&& intersection.size() * 4 > (numi + numj))
                {
                    int expected_profit = intersection.size();
                    //std::cerr << intersection.size() << " ";
                    //std::cerr << "hup #"  << (beste_parent->getSize() + intersection.size()) << "  -> "<< intersection << std::endl;
                    beste_parent->applyDownwardsInChain(this, [&](node* ancestor)
                    {
                        intersection.insert(ancestor->elements_introduced_here_relative_to_lower_chain.begin(), ancestor->elements_introduced_here_relative_to_lower_chain.end());
                        return true;
                    });
                    int before = reverse_introductions.size();

                    node* common = getNodeCreateIfMissing2(intersection, true);
                    int after = reverse_introductions.size();
                    int profit = before - after;
                    //std::cerr << before << "\t -> "  << after << "\t profit: " << profit << " expected profit: " <<  expected_profit << std::endl;
                    //assertssexec(profit >= expected_profit,  pt(beste_parent) << std::endl << pt(n) << std::endl << pt(common), dump());
                    //return n;
                }
            }

        }

        return n;
    }

    template<typename traits>
    typename traits::value_type*
    PosetOfSets<traits>::getValue(PosetNode<traits>* node)
    {
        return &node->value;
    }

    template<typename traits>
    void PosetOfSets<traits>::setValue(PosetNode<traits>* node, typename traits::value_type value)
    {
        assert(!node->isAttached());
        node->value = value;
    }

    template<typename traits>
    std::string PosetOfSets<traits>::toDot() const
    {
        //dot -Tsvg poset.dot  > poset.svg && firefox poset.svg

        std::ostringstream oss(std::ostringstream::out);
        oss << "digraph 123 { graph [splines=true,rankdir=BT]; node [shape=box]; edge [len=2, minlen=2];" << std::endl;

        std::set<node*> all_sorted;
        for(auto &i : reverse_introductions)
        {
            all_sorted.insert(i.second);
        }

        all_sorted.insert(getEmpty());
        all_sorted.insert(getAll());

        std::map<std::size_t, std::set<node*> > size_to_nodes;
        for(auto& i : all_sorted)
        {
            size_to_nodes[i->getSize()].insert(i);
        }
        for(auto& i : size_to_nodes)
        {
            if(i.second.size() > 1)
            {
                oss << "{rank=same;";
                for(auto& j : i.second)
                {
                    oss << (std::size_t) j << " ";
                }
                oss << "}\n";
            }
        }

        std::set<std::pair<size_t, size_t> > considered_edges;

        oss << (size_t) getAll() << " " << "[label=\"ALL\",style=filled,fillcolor=cyan]; { rank=sink; " << (size_t) getAll() << " }" << std::endl;
        oss << (size_t) getEmpty() << " " << "[label=\"EMPTY\",style=filled,fillcolor=\"#85FFAF\"];{ rank=source; " << (size_t) getEmpty() << " }" << std::endl;

        for(auto& n : all_sorted)
        {
            oss << (size_t) n << " " << "[label=\"" << n->toDot() << "\"];" << std::endl;
        }

        for(node* n : all_sorted)
        {
            for(node* i : n->lower)
            {
                if (!containedInVector(i->higher, n)){
                    oss << (size_t) n << "\t-> " << (size_t) i << " [constraint=false]" << std::endl;
                }
            }
            if(n->lower_chain && !containedInVector(n->lower_chain->higher, n))
            {
                oss << (size_t) n << "\t-> " << (size_t) n->lower_chain << " [color=orange, constraint=false]" << std::endl;
            }
            for(node* i : n->higher)
            {
                oss << (size_t) n << "\t-> " << (size_t) i << " [constraint=true";
                if(containedInVector(i->lower, n))
                { // bidirectional
                    oss << ",dir=both";
                }
                oss << "]" << std::endl;
            }
            for(auto& i : n->deriving_chains)
            {
                oss << (size_t) n << "\t-> " << (size_t) i << " [color=red";
                if(containedInVector(i->lower, n))
                { // bidirectional
                    oss << ",dir=both";
                }
                oss << "]" << std::endl;
            }
        }

        oss << "}" << std::endl;
        return oss.str();

    }

    template<typename traits>
    PosetNode<traits>*
    PosetOfSets<traits>::getAll() const
    {
        return (PosetNode<traits>*) &all;
    }

    template<typename traits>
    PosetNode<traits>*
    PosetOfSets<traits>::getEmpty() const
    {
        return (PosetNode<traits>*) &empty;
    }

    template<typename traits>
    PosetNode<traits>::PosetNode(const size_t size_, const PosetOfSets<traits> * myposs_, const bool attached_) :
            lower_chain(0),
            size(size_)
    {
    }

    template<typename traits>
    void PosetNode<traits>::applyDownwardsInChain(PosetOfSets<traits> * myposs, std::function<bool(PosetNode<traits>*)> f)
    {
        if(this == myposs->getEmpty())
            return;
        assert(this != myposs->getAll());

        PosetNode<traits>* h = this;
        while(h != myposs->getEmpty())
        {
            if(!f(h))
            {
                break;
            }
            h = h->lower_chain;
        }
    }

    template<typename traits>
    size_t PosetNode<traits>::getSize() const
    {
        return size;
    }

}
#endif
