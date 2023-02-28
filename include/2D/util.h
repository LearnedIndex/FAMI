#include <vector>
using namespace std;
namespace util
{
  class MyVisitor : public SpatialIndex::IVisitor
  {
  public:
    void visitNode(const SpatialIndex::INode & /* n */) override {}

    void visitData(const SpatialIndex::IData &d) override
    {
      rstar_res.push_back(d.getIdentifier());
      // the ID of this data entry is an answer to the query. I will just print it to stdout.
    }

    void visitData(std::vector<const SpatialIndex::IData *> &v) override
    {
    }
  };

  template <class KeyType>
  vector<vector<KeyType>> LoadData(const string &filename, int dim, int data_size)
  {
    vector<vector<KeyType>> data;
    ifstream in(filename);
    if (!in.is_open())
    {
      cerr << "unable to open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    vector<KeyType> item(dim);
    float key;
    while (data_size--)
    {
      if (!(in >> key))
      {
        cerr << "Not enough data entries" << endl;
        throw;
      }
      item[0] = fabs(key);
      for (int i = 1; i < dim; i++)
      {
        in >> key;
        item[i] = fabs(key);
      }
      data.push_back(item);
    }
    return data;
  }

  template <class KeyType>
  vector<vector<KeyType>> LoadData(const string &filename, int dim)
  {
    vector<vector<KeyType>> data;
    ifstream in(filename);
    if (!in.is_open())
    {
      cerr << "unable to open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    vector<KeyType> item(dim);
    float key;
    while (in >> key)
    {
      item[0] = fabs(key);
      for (int i = 1; i < dim; i++)
      {
        in >> key;
        item[i] = fabs(key);
      }
      data.push_back(item);
    }
    return data;
  }

  template <class KeyType>
  vector<KeyType> GetMinKey(const vector<vector<KeyType>> &data, int dim)
  {
    vector<KeyType> min_keys(dim);
    int size = data.size();
    for (int i = 0; i < dim; i++)
    {
      min_keys[i] = data[0][i];
      for (int j = 1; j < size; j++)
      {
        min_keys[i] = min(min_keys[i], data[j][i]);
      }
    }
    return min_keys;
  }
  template <class KeyType>
  vector<KeyType> GetMaxKey(const vector<vector<KeyType>> &data, int dim)
  {
    vector<KeyType> max_keys(dim);
    int size = data.size();
    for (int i = 0; i < dim; i++)
    {
      max_keys[i] = data[0][i];
      for (int j = 1; j < size; j++)
      {
        max_keys[i] = max(max_keys[i], data[j][i]);
      }
    }
    return max_keys;
  }

  float gridsize[2];
  template <class KeyType>
  vector<vector<KeyType>> GetUpperData(const vector<vector<KeyType>> &data, int dim, vector<float> max_keys)
  {
    std::default_random_engine e;
    std::normal_distribution<float> n(0, 1);
    vector<vector<float>> upper_data;
    vector<KeyType> item(dim);
    int size = data.size();
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < dim; j++)
      {
        // item[j] = data[i][j] + max(0.0, gridsize[j] + n(e));
        item[j] = min(data[i][j] + gridsize[j], max_keys[j]);
      }
      upper_data.push_back(item);
    }
    return upper_data;
  }

  float gridsizeX;
  template <class KeyType>
  vector<KeyType> GetUpperDataX(const vector<vector<KeyType>> &data, int dim, vector<float> max_keys)
  {
    vector<float> upper_dataX;
    int size = data.size();
    for (int i = 0; i < size; i++)
    {
      upper_dataX.push_back(min(data[i][0] + gridsizeX, max_keys[0]));
    }
    return upper_dataX;
  }

  float gridsizeY;
  template <class KeyType>
  vector<KeyType> GetUpperDataY(const vector<vector<KeyType>> &data, int dim, vector<float> max_keys)
  {
    vector<float> upper_dataY;
    int size = data.size();
    for (int i = 0; i < size; i++)
    {
      upper_dataY.push_back(min(data[i][1] + gridsizeY, max_keys[1]));
    }
    return upper_dataY;
  }
}