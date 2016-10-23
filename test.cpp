#include "posetofsets.h"

struct PNTraits
{
    typedef std::size_t subelement_type;
    typedef std::size_t value_type;
    typedef std::set<std::size_t> subelementlookup;
//    typedef subelementlookup rarestfirst;
};

void test_empty_search()
{
    PosetOfSets<PNTraits> poss;
    typedef PosetNode<PNTraits> node;
    PNTraits::subelementlookup mmm;
    PosetQuery<PNTraits> query(&poss);
    std::set<node*> largest_subsets, smallest_sets_surrounding;
    query.find_largest_subsets_of_given_set(mmm, largest_subsets);

    assert2(largest_subsets.size() == 1, largest_subsets.size());
    assert3(*largest_subsets.begin() == poss.getEmpty(), largest_subsets, poss.getEmpty());
    query.find_smallest_sets_surrounding_given_set(mmm, smallest_sets_surrounding);

    assert3(smallest_sets_surrounding == std::set<node*>
    { poss.getEmpty() }, smallest_sets_surrounding, poss.getEmpty());

    mmm.insert(5);
    smallest_sets_surrounding.clear();
    query.find_smallest_sets_surrounding_given_set(mmm, smallest_sets_surrounding);
    assert3(smallest_sets_surrounding == std::set<node*>
    { poss.getAll() }, smallest_sets_surrounding, poss.getAll());

    node* five = poss.getNodeCreateIfMissing2(mmm);
    assert2(five && five != poss.getAll() && five != poss.getEmpty(), five);
    assert2(five->getSize() == 1, five->getSize());
    assertssexec(five->elements_introduced_here_relative_to_lower_chain.size() == 1, pt(five->elements_introduced_here_relative_to_lower_chain.size()), writeToFile("/tmp/badposet.dot", poss.toDot()));

    assertssexec(poss.getEmpty()->deriving_chains == std::vector<node*>
    { five }, pt(poss.getEmpty()->deriving_chains), poss.dump());
    assert(poss.getEmpty()->higher == std::vector<node*>
    { five });
    assert(poss.getAll()->lower_chain == 0);
    assert(poss.getAll()->lower == std::vector<node*>
    { five });

    assert(poss.reverse_introductions.size() == 1);
    assert(poss.reverse_introductions.find(5) != poss.reverse_introductions.end());
    assert(poss.reverse_introductions.find(5)->second == five);

    smallest_sets_surrounding.clear();
    query = PosetQuery<PNTraits>(&poss);
    query.find_smallest_sets_surrounding_given_set(mmm, smallest_sets_surrounding);

    assertss(smallest_sets_surrounding == std::set<node*>
    { five }, smallest_sets_surrounding << " != " << std::set<node*> {five} << " empty: " << poss.getEmpty() << " all: " << poss.getAll());

    largest_subsets.clear();
    query.find_largest_subsets_of_given_set(mmm, largest_subsets);
    assert2(largest_subsets.size() == 1, largest_subsets.size());

    assertss(largest_subsets == std::set<node*>
    { five }, largest_subsets << " != " << five << " empty: " << poss.getEmpty() << " all: " << poss.getAll());

    node* five_again = poss.getNodeCreateIfMissing2(mmm);
    assert2(five_again && five_again != poss.getAll() && five_again != poss.getEmpty(), five_again);
    assert2(five_again->getSize() == 1, five_again->getSize());
    assertss(five == five_again, pt(five_again));

    mmm =
    {   1,2,3};

    node* onetwothree = poss.getNodeCreateIfMissing2(mmm);
    assert(onetwothree);
    assert(onetwothree->getSize() == 3);
    assert(onetwothree->elements_introduced_here_relative_to_lower_chain == mmm);

    mmm =
    {   1,2,3,4,5};

    smallest_sets_surrounding.clear();
    query = PosetQuery<PNTraits>(&poss);
    query.find_smallest_sets_surrounding_given_set(mmm, smallest_sets_surrounding);
    assertss(smallest_sets_surrounding.size() == 1, pt(smallest_sets_surrounding) << pt(poss.getAll()) << pt(poss.getEmpty()) << pt(five) << pt(onetwothree));
    assertss(smallest_sets_surrounding == std::set<node*>
    { poss.getAll() }, smallest_sets_surrounding << " != " << poss.getAll() << " empty: " << poss.getEmpty() << " all: " << poss.getAll());

    node* onetwothreefourfive = poss.getNodeCreateIfMissing2(mmm);
    assert(onetwothreefourfive);
    assert2(onetwothreefourfive->getSize() == 5, onetwothreefourfive->toString());

    assert(onetwothreefourfive->elements_introduced_here_relative_to_lower_chain.size() == 2);
}

void writeToFile(const std::string file, const std::string& s)
{
    FILE * f = fopen(file.c_str(), "w");
    if(f != NULL)
    {
        fputs(s.c_str(), f);
        fclose(f);
    }
}

// probleem: stel ik maak echt uitsluitend nog de links die duidelijk worden door boven en onder-sets te berekenen ...
// kan ik dan nog gemakkelijk die chains houden?
// stel ik wil die regel houden dat een chain die erft van chain C ook later moet convergeren in chain x ...
// ... ik ga moeten toelaten dat het divergeren niet letterlijk chain C mag zijn maar ook een voorouder ...
// ... en bij het convergeren ... een divergente tak ... van waar het van komt ...
// die chaintoestanden lijken vooral nodig om ... te weten waar de introduced elements bij komen ...
//
//
// en de 2 algos ... wat hebben die precies nodig .. qua constraints ...
// is het relevant om te weten wat er convergeert?? laat de chains gewoon doodlopen naar boven toe ? (de higher blijven bestaan maar die higher_chain moet dan uitstaan somehow)
// Ah, om te voorkomen dat die introduced-sets veel te groot worden, die kunnen veel kleiner door te convergeren ...
// dus ... als de higher_chain != my_chain ... dan euhm ... indiceert dat convergentie ... hoe hoort dan die introduced-chain te werken ???? gaat toch niet??
//
//
// en hoe verloopt die addition van een nieuwe chain ...
//
// die chain-index ... wordt dan toch slechts een aftakking-index...
//

// Stel:
//  123457   12346
//    |   /
//   123
//
// waar hoort   1234? ?
// tussen 123 en 123457  lijkt voor de hand liggend wegens zelfde chainindex ...
// de keuze ...
//  voor de onderliggende node best die node pakken met de grootste set ... zodat de introduced set klein kan gehouden worden...
// voor wat daarboven betreft ... in dit voorbeeld moet die divergentie naar boven gepusht worden ...
// het doel is dan ook ... om zoveel mogelijk van dergelijke divergenties ... omhoog te duwen ...
//  want voor elk zulks geval moet er minder data redundant geintroduceerd worden ...
//
// ... die bovenset moet gepartitioneerd worden volgens de elementen in de onderset!
// dus: A. zo hoog mogelijke node onder de nieuwe ... en B. zoveel mogelijk nodes laten divergeren vanaf de nieuwe door ze te pushen naar boven...
// ... maar .. uit die hogere set ... tellen precies alleen ... die direct met de lagere verbonden zijn ...
//
// en als er nergens een directe verbinding is tussen onder en boven om tussen te duiken ...
// dan moet ... er een branch gemaakt worden ... ...
// die totaal niet moet convergeren ...
// en best vanaf een zo hoog mogelijke node ... kwestie zo weinig mogelijk elementen te moeten introduceren ...
// voila ... voila ...

//void test_huge_posets

void create_nice_graph(std::size_t max_elements, std::size_t max_nodes)
{
    PosetOfSets<PNTraits> poss;
    std::size_t total_subelements = 0;
    std::vector<std::set<std::size_t> > contents;
    for(std::size_t i = 0; i < max_nodes; ++i)
    {
        std::set<std::size_t> contents;

        std::size_t numelements = rand() % max_elements;

        while(contents.size() < numelements)
        {
            std::size_t base = (rand() % max_elements);
            //for (int j = 0; j < 10; ++j){
            contents.insert(base  );
            //}
        }
        total_subelements += contents.size();
        poss.getNodeCreateIfMissing2(contents);
        if (i % 1000 == 0){
            prt4(max_elements, poss.reverse_introductions.size(), max_nodes, total_subelements);
        }
    }


   poss.dump();

}

void test_random_poset_nodes(std::size_t max_elements, std::size_t max_nodes)
{
    //std::size_t max_elements = 6;
    //std::size_t max_nodes = 2000;
    std::size_t total_subelements = 0;
    std::vector<std::set<std::size_t> > contents;

    for(std::size_t i = 0; i < max_nodes; ++i)
    {
        contents.push_back(std::set<std::size_t>());
        std::size_t numelements = rand() % max_elements;
//        for(std::size_t j = 0; j < numelements; ++j)
        while(contents.rbegin()->size() < numelements)
        {
            contents.rbegin()->insert((rand() % max_elements) + 100);
        }
        contents.rbegin()->insert(1000);
        contents.rbegin()->insert(1001);
        contents.rbegin()->insert(1003);
        total_subelements += contents.rbegin()->size();
    }

    PosetOfSets<PNTraits> poss;
    typedef PosetNode<PNTraits> node;
    std::vector<node*> nodes;
    std::set<node*> nodes_sorted;
    //std::map<std::size_t, std::map<std::size_t, node*> > nodegrid;

    for(std::size_t i = 0; i < max_nodes; ++i)
    {

        // prt2(i, contents[i]);
        //std::cerr << "[" << contents[i].size();
        node* c = poss.getNodeCreateIfMissing2(contents[i]);
        nodes.push_back(c);
        nodes_sorted.insert(c); // to count unique nodes.
        //std::cerr << "]";
        if(true)
        {
//        writeToFile("/tmp/badposet_" + std::to_string(i) + ".dot", poss.toDot());

            //assertssexec(nodegrid[c->getChain()][c->getSize()] == 0 || nodegrid[c->getChain()][c->getSize()] == c, pt(nodegrid[c->getChain()][c->getSize()]) << std::endl << pt(c) << std::endl <<  pt(c->toDot()) << std::endl << contents[i], poss.dump());
//        nodegrid[c->getChain()][c->getSize()] = c;

            for(auto& higher : c->higher)
            {
                assert(higher->getSize() > c->getSize());
            }

            for(auto& lower : c->lower)
            {
                assert(lower->getSize() < c->getSize());
            }

//        for(auto &j : nodegrid)
//        {
//            for(auto &k : j.second)
//            {
//                assertssexec(!k.second->elements_introduced_here_relative_to_lower_chain.empty(), pt(contents[i] ), poss.dump());
//            }
//        }

            assertssexec(contents[i].empty() || !c->elements_introduced_here_relative_to_lower_chain.empty(), pt(contents[i]), poss.dump());

            //writeToFile("/tmp/badposetgrr.dot", poss.toDot());

//        assertssexec(
//                nodegrid[c->getChain()].begin()->second->getSize() == nodegrid[c->getChain()].begin()->second->elements_introduced_here_relative_to_lower_chain.size()
//                        || nodegrid[c->getChain()].begin()->second->lower_chain != c->getChain(),
//                (nodegrid[c->getChain()].begin()->second->getSize() == nodegrid[c->getChain()].begin()->second->elements_introduced_here_relative_to_lower_chain.size()) << std::endl << pt(c->getChain()) << std::endl << pt(c->toDot()) << std::endl << pt(contents[i]),
//                poss.dump());

//        assertssexec(nodegrid[c->getChain()].rbegin()->second->chain_continues == false, pt(nodegrid[c->getChain()].rbegin()->second->toDot()), poss.dump());

//        for(auto j = nodegrid[c->getChain()].begin(); j != nodegrid[c->getChain()].end(); ++j)
//        {
//            assert(j->first == nodegrid[c->getChain()].rbegin()->first || j->second->chain_continues == true);
//            assertssexec(j == nodegrid[c->getChain()].begin() || j->second->lower_chain == c->getChain(), "blaat", poss.dump());
//        }
            // als er hoger of lager een element is in die chain ... dan moet dat opvallen in hun continue_chain en lower_chain

            std::map<std::size_t, std::map<std::size_t, node*> > element_to_chains_to_node;
            for(auto& j : poss.reverse_introductions)
            {
//            assertssexec(element_to_chains_to_node[j.first][j.second->getChain()] == 0, pt(poss.reverse_introductions), poss.dump());
                //           element_to_chains_to_node[j.first][j.second->getChain()] = j.second;
                assert(j.second->elements_introduced_here_relative_to_lower_chain.find(j.first) != j.second->elements_introduced_here_relative_to_lower_chain.end());
            }

            assert(c->getSize() == contents[i].size());

            for(std::size_t j = 0; j <= i; ++j)
            {
                for(auto& subsets : nodes[j]->lower)
                {
                    for(node* one_above_that_highest_subset : subsets->higher)
                    {
                        if(one_above_that_highest_subset->lower_chain == subsets && one_above_that_highest_subset != poss.getAll() && one_above_that_highest_subset != nodes[j])
                        {
                            assertssexec(
                                    !std::includes(contents[j].begin(), contents[j].end(), one_above_that_highest_subset->elements_introduced_here_relative_to_lower_chain.begin(),
                                            one_above_that_highest_subset->elements_introduced_here_relative_to_lower_chain.end()),
                                    pt(max_elements) << std::endl << pt(i) << std::endl << pt(one_above_that_highest_subset) << std::endl << pt(nodes[j]), poss.dump());
                        }
                    }
                }
            }

            std::set<std::size_t> reconstructed;
            c->applyDownwardsInChain(&poss, [&](node* ancestor)
            { //prt5(ancestor, ancestor->getSize(), ancestor->lower, ancestor->lower_chain, poss.getEmpty()); //
                        reconstructed.insert(ancestor->elements_introduced_here_relative_to_lower_chain.begin(),ancestor->elements_introduced_here_relative_to_lower_chain.end()); return true;});

            assertssexec(reconstructed == contents[i], pt(nodes[i]) << std::endl << pt(reconstructed) << std::endl << pt(contents[i]) << std::endl << pt(poss.reverse_introductions), poss.dump());

            //      assertssexec(i != 20, pt(reconstructed) << std::endl << pt(contents[i])     , writeToFile("/tmp/badposet.dot", poss.toDot()));
        }
    }

    for(std::size_t i = 0; i < max_nodes; ++i)
    {
//        prt2(i, contents[i]);
        node* c = poss.getNodeCreateIfMissing2(contents[i]);
        assert(c->getSize() == contents[i].size());

        std::set<std::size_t> reconstructed;
        c->applyDownwardsInChain(&poss, [&](node* ancestor)
        {   reconstructed.insert(ancestor->elements_introduced_here_relative_to_lower_chain.begin(),ancestor->elements_introduced_here_relative_to_lower_chain.end()); return true;});
        assertssexec(reconstructed == contents[i], pt(reconstructed) << std::endl << pt(contents[i]), poss.dump());
        assertssexec(c == nodes[i], c << std::endl << nodes[i], poss.dump());

    }

    prt7(max_elements, poss.reverse_introductions.size(), nodes_sorted.size(), max_nodes, total_subelements, contents[0].size(), contents[1].size());
//    assertssexec(poss.reverse_introductions.size() > 3,pt(poss.reverse_introductions),poss.dump());
//    writeToFile("/tmp/fun_"".dot", poss.toDot());
    //assertssexec(false, "", poss.dump());
}

void test_common_node()
{
    PosetOfSets<PNTraits> poss;
    std::set<std::size_t> a
    { 1, 2, 3, 4, 5, 6 };
    std::set<std::size_t> b
    { 1, 2, 3, 4, 5, 7 };
    poss.getNodeCreateIfMissing2(a, true);
    poss.getNodeCreateIfMissing2(b, true);
    poss.dump();
    assert(poss.reverse_introductions.size() == 7);
}

int main()
{
    srand(4);
    create_nice_graph(4,10);
    abort();

    test_empty_search();
    //test_common_node();

    //test_random_poset_nodes(100, 20);

    //test_random_poset_nodes(6,200);
    for(int i = 2; i < 128; ++i)
    {
        for(int seed = 0; seed < 1; ++seed)
        {
            //prt2(i, seed);
            //int seed = 48;
            srand(seed);
            //prt(seed); //pow(2,i)
            prt2(i,seed);
            test_random_poset_nodes(pow(2,i),100);
        }
    }
    std::cerr << std::endl << "done" << std::endl;

//    mmm.insert(5);
//    poss.getNodeCreateIfMissing2(mmm);

}
