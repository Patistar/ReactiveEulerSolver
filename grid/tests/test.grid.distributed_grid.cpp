#include "hpx/hpx.hpp"

class partition_data {
public:
  explicit partition_data(std::size_t size);

private:
  std::vector<double> m_data;

  // Serialization support: even if all of the code below runs on one
  // locality only, we need to provide an (empty) implementation for the
  // serialization as all arguments passed to actions have to support this.
  friend class hpx::serialization::access;

  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version) {
    (void)(ar & m_data);
  }
};

partition_data::partition_data(std::size_t size)
    : m_data(std::vector<double>(size, 0.0)) {}

int hpx_main() { return hpx::finalize(); }

int main(int argc, char** argv) { return hpx::init(argc, argv); }
