#ifndef NODE
#define NODE

#include <immintrin.h>

/*! \brief Base class defines searching interfaces
 */
class Node {
	public:
		Node(std::vector<float> feature_vector,
				uint32_t id,
				uint32_t delimiter_dimension) :
			data_object(feature_vector),
			id(id),
			delimiter_dimension(delimiter_dimension) {
				left_child = NULL;
				right_child = NULL;
			};
		virtual ~Node() {};

		int insertObject(std::vector<float> feature_vector, uint32_t id, int depth) {
                        depth++;
			if (feature_vector[delimiter_dimension] <= data_object[delimiter_dimension]) {
				if (left_child == NULL) {
					left_child = new Node(feature_vector, id, (delimiter_dimension + 1) % data_object.size());
				} else {
					depth = left_child->insertObject(feature_vector, id, depth);
				}
			} else {
				if (right_child == NULL) {
					right_child = new Node(feature_vector, id, (delimiter_dimension + 1) % data_object.size());
				} else {
					depth = right_child->insertObject(feature_vector, id, depth);
				}
			}

			return depth;
		}

		std::vector<uint32_t> rangeSearch(std::vector<uint32_t> &results, std::vector<float> min, std::vector<float> max) {
			bool match = true;
			for (size_t i = 0; i < data_object.size(); ++i) {
				match = (match && data_object[i] >= min[i] && data_object[i] <= max[i]);
			}
			if (match)
				results.push_back(id);
	
			if (left_child != NULL && min[delimiter_dimension] <= data_object[delimiter_dimension])
				left_child->rangeSearch(results, min, max);	
		 	if(right_child != NULL && max[delimiter_dimension] > data_object[delimiter_dimension])
				right_child->rangeSearch(results, min, max);	

			return results;
		}

		std::vector<uint32_t> rangeSearchSIMD(std::vector<uint32_t> &results, std::vector<float> min, std::vector<float> max) {
			bool match = true;
			size_t i = 0;
			__m256 lower_reg, upper_reg, search_reg, lower_res, upper_res;
			size_t mask_lower, mask_upper, mask;
			size_t simd_compares = ((size_t) (data_object.size() / 8)) * 8;
			for (; i < simd_compares; i += 8) {
				lower_reg = _mm256_loadu_ps(&min[i]);
				upper_reg = _mm256_loadu_ps(&max[i]);
				search_reg = _mm256_loadu_ps(&this->data_object[i]);
				lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
				mask_lower = _mm256_movemask_ps(lower_res);
				upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
				mask_upper = _mm256_movemask_ps(upper_res);
				mask       = mask_lower & mask_upper;

				// mask is set to 0xFF if every SIMD lane matchs
				if (mask < 0xFF) {
					match = false;
					i = data_object.size();
					break;
				}
			}
			for (; i < data_object.size(); ++i) {
				if (data_object[i] < min[i] || data_object[i] > max[i]) {
					match = false; break;
				}
			}
			if (match)
				results.push_back(id);

			if (left_child != NULL && min[delimiter_dimension] <= data_object[delimiter_dimension])
				left_child->rangeSearch(results, min, max);	
			if(right_child != NULL && max[delimiter_dimension] > data_object[delimiter_dimension])
				right_child->rangeSearch(results, min, max);	

			return results;
		}

	private:
		std::vector<float> data_object;
		uint32_t id;
		uint32_t delimiter_dimension;
		Node* left_child;
		Node* right_child;
};
#endif // NODE
