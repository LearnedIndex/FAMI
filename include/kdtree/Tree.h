#ifndef KDTREE
#define KDTREE

#include "Node.h"

#include <vector>

class Tree {
public:

    Tree() {
      root_node = NULL;
      depth = 0;
    }

    ~Tree() {
      delete root_node;
    }

    void insertObject(std::vector<float> feature_vector, uint32_t id) {
	    if (root_node == NULL) {
		    root_node = new Node(feature_vector, id, 0);
	    } else {
		    int new_depth = root_node->insertObject(feature_vector, id, 0);
		    if (this->depth < new_depth)
		      this->depth = new_depth;
	    }
    }

    std::vector<uint32_t> rangeSearch(std::vector<float> min, std::vector<float> max) {
	    std::vector<uint32_t> results;

	    if (root_node != NULL)
		    results = root_node->rangeSearch(results, min, max);

	    return results;
    }

    std::vector<uint32_t> rangeSearchSIMD(std::vector<float> min, std::vector<float> max) {
	    std::vector<uint32_t> results;

	    if (root_node != NULL)
		    results = root_node->rangeSearchSIMD(results, min, max);

	    return results;
    }

    int getDepth() {
    	return this->depth;
    }

private:
    Node *root_node;
    int depth;
};

#endif // KDTREE
