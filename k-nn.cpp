/*
 * This reads in a binary file of doubles, and inserts into a BST.
 */

#include <random>
#include <cstdlib>
#include <cstdio>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <sys/mman.h>
#include <linux/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <ostream>
#include <math.h>
#include <thread>
uint64_t points;
int dimensions;
int q_points;
uint64_t train_id;
uint64_t n_threads;

class Node {
    public:
        Node(std::vector<float> v,int l) : value{v},level{l}, left_child{nullptr}, right_child{nullptr} {}
        static void insert(Node **pp, Node *n);
        static void verify(Node *n);
        static void knn_algo(std::vector<float> &,uint64_t ,Node **p);
        static float distance_calculator(std::vector<float> &,std::vector<float> &);
    private:
        std::vector<float> value;
        int level;
        Node *left_child, *right_child;
};
Node *root = nullptr;
Node *nodes = nullptr;
struct Near{
  std::vector<float> cluster;
  float distance;
};

void
Node::insert(Node **pp, Node *n) {
    if (*pp == nullptr) {
        *pp = n;
                // printf("tree:%f\n",n->value[n->level]);
    } else {
        Node *p = *pp;

        if (n->value[p->level] <= p->value[p->level]) {
            insert(&p->left_child,n);
        } else {
            insert(&p->right_child,n);
        }
    }
    return;
}

void
Node::verify(Node *n) {
    if (n->left_child) {
        assert(n->left_child->value[n->level] <= n->value[n->level]);
        verify(n->left_child);
    }
    for(int i=0;i<dimensions;i++)
    {
    // printf("       %f             %d",n->value[i],n->level);
    }
    // printf("\n");
    if (n->right_child) {
        assert(n->right_child->value[n->level] > n->value[n->level]);
        verify(n->right_child);
    }
    return;
}

class Reader {
    public:
        Reader(const char *p) : ptr{p} {}
        template <typename T>
        Reader &operator>>(T &o) {
            // Assert alignment.
            assert(uintptr_t(ptr)%sizeof(T) == 0);
            o = *(T *) ptr;
            ptr += sizeof(T);
            return *this;
        }
    private:
        const char *ptr;
};

void dump(const std::string &fn);

void split(std::vector<std::vector<float> > &, uint64_t ,uint64_t ,int ,uint64_t );

std::vector<std::vector<float> > val;
std::vector<int> split_level;
void split(std::vector<std::vector<float> > &array,uint64_t l,uint64_t r,uint64_t level,uint64_t d)
{
    if(array.size()!=0)
    {
    if(level==d)
    {
      level =0;
    }
    std::sort(array.begin(),array.end(), [=](const std::vector<float> &p1,const std::vector<float> &p2){return (p1[level] < p2[level]);});
    val.push_back(array.at(array.size()/2));
    split_level.push_back(level);
    array.erase(array.begin() + array.size()/2);
    level++;
    l=0;
    uint64_t m = array.size()/2;
    r =array.size();
	  std::vector<std::vector<float> > l_array;
    std::vector<std::vector<float> > r_array;
              std::vector<float> row;
    for (uint64_t i = l; i < m; i++) {
            row =array.at(i);
        l_array.push_back(row);
    }
    for (uint64_t i= m; i < r; i++) {

            row =array.at(i);
        r_array.push_back(row);
    }

    split(l_array, l, m,level,d);
    split(r_array, m, r,level,d);
    }
}


float Node::distance_calculator(std::vector<float> &point1,std::vector<float> &point2)
{
  float distance=0;
  for(uint64_t i=0; i<point1.size();i++)
  {
    distance=distance+pow((point1[i]-point2[i]),2);
  }
  distance =sqrt(distance);
  // printf("Distance:%f\n",distance);
  return distance;
}
std::vector<Near> n_points;
float c_dist;
float l_dist;
void Node::knn_algo(std::vector<float> &point,uint64_t k,Node **pp) {
  if (*pp == nullptr) {
      return;
  }
  Node *p = *pp;
  c_dist = Node::distance_calculator(point,p->value);
  std::sort(n_points.begin(),n_points.end(), [=](const Near &p1,const Near &p2){return (p1.distance < p2.distance);});
  if(n_points.size() !=0)
  {
    l_dist = n_points[n_points.size()-1].distance;
  }
  if(n_points.size()<k)
  {
    n_points.push_back({p->value,c_dist});
  }
  else if(l_dist>c_dist)
  {
    n_points.pop_back();
    n_points.push_back({p->value,c_dist});
  }
  // else if(l_dist<=c_dist)
  // {
  //   return;
  // }
  Node *near_tree;
  Node *far_tree;
      if (point[p->level] <= p->value[p->level])
      {
        near_tree = p->left_child;
      }
      else
      {
        near_tree = p->right_child;
      }
      if (point[p->level] > p->value[p->level])
      {
        far_tree = p->left_child;
      }
      else
      {
        far_tree = p->right_child;
      }

      knn_algo(point,k,&near_tree);
      knn_algo(point,k,&far_tree);
}


int
main(int argc, char **argv) {

  if (argc != 5) {
      std::fprintf(stderr, "./k-nn <n_cores> <training_file> <query_file> <result_file>\n");
      std::exit(1);
      n_threads = (uint64_t)argv[1];
  }
    for (int i = 2; i < argc; i++) {
        dump(argv[i]);
    }
}

void
dump(const std::string &fn) {

    std::cout << fn << std::endl;

    /*
     * Use mmap() for convenience.
     */

    int fd = open(fn.c_str(), O_RDONLY);
    if (fd < 0) {
        int en = errno;
        std::cerr << "Couldn't open " << fn << ": " << strerror(en) << "." << std::endl;
        exit(2);
    }

    // Get the actual size of the file.
    struct stat sb;
    int rv = fstat(fd, &sb); assert(rv == 0);
    // std::cout << sb.st_size << std::endl;

    // Use some flags that will hopefully improve performance.
    void *vp = mmap(nullptr, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (vp == MAP_FAILED) {
        int en = errno;
        fprintf(stderr, "mmap() failed: %s\n", strerror(en));
        exit(3);
    }
    char *file_mem = (char *) vp;

    // Tell the kernel that it should evict the pages as soon as possible.
    rv = madvise(vp, sb.st_size, MADV_SEQUENTIAL|MADV_WILLNEED); assert(rv == 0);

    rv = close(fd); assert(rv == 0);

    // Prefix to print before every line, to improve readability.
    std::string pref("    ");

    /*
     * Read file type string.
     */
    int n = strnlen(file_mem, 8);
    std::string file_type(file_mem, n);
    std::cout << pref << "File type string: " << file_type << std::endl;

    // Start to read data, skip the file type string.
    Reader reader{file_mem + 8};

    // TODO: Code below is repetitive, cleanup.
    // TODO: Add io manip to print with commas.

    if (file_type == "TRAINING") {

        uint64_t id;
        uint64_t n_points;
        uint64_t n_dims;
        reader >> id >> n_points >> n_dims;
        points = n_points;
        dimensions = n_dims;
        train_id = id;
        printf("Building tree:-\n");
        std::vector<std::vector<float> > a_points;
        for (std::uint64_t i = 0; i < n_points; i++) {
              std::vector<float> row;
            for (std::uint64_t j = 0; j < n_dims; j++) {
                float f;
                reader >> f;
                row.push_back(f);
            }
            a_points.push_back(row);
        }
        // for (int i = 0; i < a_points.size(); i++)
        // {
        //   for (int j = 0; j < a_points[i].size(); j++)
        //   {
        //     // std::cout << a_points[i][j];
        //   // printf("          ");
        //   }
        //
        //           // printf("\n");
        // }
                          // printf("\n");
        split(a_points,0,n_points,(uint64_t)0,n_dims);
        // print("Hello");
        nodes = (Node *) ::operator new(points*sizeof(Node));
        for (uint64_t i = 0; i < n_points; i++) {
          Node *n = new (nodes +i) Node{val.at(i),split_level.at(i)};
          Node::insert(&root, n);
        }
    } else if (file_type == "QUERY") {

        uint64_t id;
        uint64_t n_queries;
        uint64_t n_dims;
        uint64_t n_neighbors;
        printf("Searching Neighbouring points for each query:-\n");
        reader >> id >> n_queries >> n_dims >> n_neighbors;
        q_points=n_neighbors;
        std::cout << pref << "Query file ID: " << std::hex << std::setw(16) << std::setfill('0') << id << std::dec << std::endl;
        std::cout << pref << "Number of queries: " << n_queries << std::endl;
        std::cout << pref << "Number of dimensions: " << n_dims << std::endl;
        std::cout << pref << "Number of neighbors to return for each point: " << n_neighbors << std::endl;
        std::vector<std::vector<float> > query_points;
        for (std::uint64_t i = 0; i < n_queries; i++) {
            // std::cout << pref << "Query " << i << ": ";
            std::vector<float> row;
            for (std::uint64_t j = 0; j < n_dims; j++) {
                float f;
                reader >> f;
                row.push_back(f);
            }
            query_points.push_back(row);
            // std::cout << std::endl;
        }

        // for (int i = 0; i < query_points.size(); i++)
        // {
        //   for (int j = 0; j < query_points[i].size(); j++)
        //   {
        //     std::cout << query_points[i][j];
        //   // printf("          ");
        //   }
        //
        //           // printf("\n");
        // }
                          // printf("\n");
        std::ofstream file;
        file.open("result_file", std::ios::binary);
        if (!file)
        {
          int en = errno;
          std::fprintf(stderr, "Couldn't open %s: %s\n", "result_file", strerror(en));
          exit(4);
        }
        uint64_t result_id =rand()%10000;
        file.write("RESULT", sizeof(char)*8);
        file.write(reinterpret_cast<char *>(&train_id), sizeof(uint64_t));
        file.write(reinterpret_cast<char *>(&id), sizeof(uint64_t));
        file.write(reinterpret_cast<char *>(&result_id), sizeof(uint64_t));
        file.write(reinterpret_cast<char *>(&n_queries), sizeof(uint64_t));
        file.write(reinterpret_cast<char *>(&n_dims), sizeof(uint64_t));
        file.write(reinterpret_cast<char *>(&n_neighbors), sizeof(uint64_t));
        for(uint64_t i=0;i<n_queries;i++)
        {
            Node::knn_algo(query_points[i],q_points,&root);
                        // printf("Query %d nearest points\n",i);
                        // for(int i=0;i<n_points.size();i++)
                        // {
                        //   for (int k = 0; k < n_dims; k++)
                        //     {
                        //   // printf("X dim:%f and its distance:%f \n",n_points[i].cluster[k],n_points[i].distance);
                        // }
                        // }
            for (uint64_t j = 0; j < n_points.size(); j++)
            {
                for (uint64_t k = 0; k < n_dims; k++)
                {
                      float x = n_points[j].cluster[k];
                      file.write(reinterpret_cast<char *>(&x), sizeof(x));
                                  // printf("Val:%f\n",x);
                }
            }
                            // std::cout << x << std::endl;
            n_points.clear();
        }
      	file.close();
        // std::cerr << "Verifying..." << std::endl;
        Node::verify(root);
        for (uint64_t i = 0; i < points; i++) {
            nodes[i].~Node();
        }

        ::operator delete(nodes);
    } else if (file_type == "RESULT") {

        uint64_t training_id;
        uint64_t query_id;
        uint64_t result_id;
        uint64_t n_queries;
        uint64_t n_dims;
        uint64_t n_neighbors;

        reader >> training_id >> query_id >> result_id >> n_queries >> n_dims >> n_neighbors;
        printf("To print result please remove comments on print\n");
        std::cout << pref << "Training file ID: " << std::hex << std::setw(16) << std::setfill('0') << training_id << std::dec << std::endl;
        std::cout << pref << "Query file ID: " << std::hex << std::setw(16) << std::setfill('0') << query_id << std::dec << std::endl;
        std::cout << pref << "Result file ID: " << std::hex << std::setw(16) << std::setfill('0') << result_id << std::dec << std::endl;
        std::cout << pref << "Number of queries: " << n_queries << std::endl;
        std::cout << pref << "Number of dimensions: " << n_dims << std::endl;
        std::cout << pref << "Number of neighbors returned for each query: " << n_neighbors << std::endl;
        // for (std::uint64_t i = 0; i < n_queries*n_neighbors; i++) {
        //     std::cout << pref << "Result " << i << ": ";
        //     for (std::uint64_t j = 0; j < n_dims; j++) {
        //         float f;
        //         reader >> f;
        //         std::cout << std::fixed << std::setprecision(6) << std::setw(15) << std::setfill(' ') << f;
        //         // Add comma.
        //         if (j < n_dims - 1) {
        //             std::cout << ", ";
        //         }
        //     }
        //     std::cout << std::endl;
        // }

    } else {
        std::cerr << "Unknown file type: " << file_type << std::endl;
        exit(2);
    }

    rv = munmap(file_mem, sb.st_size); assert(rv == 0);

    // For this simple struct, it's not really necessary to call the dtor, but
    // just some defensive programming.

    struct rusage ru;
    rv = getrusage(RUSAGE_SELF, &ru); assert(rv == 0);
    auto cv = [](const timeval &tv) {
        return double(tv.tv_sec) + double(tv.tv_usec)/1000000;
    };
    std::cerr << "Resource Usage:\n";
    std::cerr << "    User CPU Time: " << cv(ru.ru_utime) << '\n';
    std::cerr << "    Sys CPU Time: " << cv(ru.ru_stime) << '\n';
    std::cerr << "    Max Resident: " << ru.ru_maxrss << '\n';
    std::cerr << "    Page Faults: " << ru.ru_majflt << '\n';



}
