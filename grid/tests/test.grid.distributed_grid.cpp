#include <hpx/hpx.hpp>

class component_server
    : public hpx::components::component_base<component_server> {
public:
  component_server() = default;
};

using component_name = hpx::components::component<component_server>;
HPX_REGISTER_COMPONENT(component_name);

class component_client
    : public hpx::components::client_base<component_client, component_server> {
  using base_type =
      hpx::components::client_base<component_client, component_server>;
public:
  component_client() = default;
  component_client(hpx::id_type where)
      : base_type(hpx::new_<component_server>(std::move(where))) {}
  component_client(hpx::future<component_client> node)
      : base_type(std::move(node)) {}
  component_client(hpx::future<hpx::id_type>&& id) : base_type(std::move(id)) {}
};

struct function_type {
  static void invoke() {}
};

struct action_type
    : hpx::actions::make_action<decltype(&function_type::invoke),
                                &function_type::invoke, action_type> {};

int hpx_main() {
  const std::vector<hpx::id_type> localities = hpx::find_all_localities();
  component_client node(localities[0]);
#ifdef NO_DATAFLOW
  hpx::async(action_type(), node.get_id()).get();
#else
  hpx::dataflow(action_type(), node.get_id()).get();
#endif
  return hpx::finalize();
}

int main(int argc, char** argv) { hpx::init(argc, argv); }
