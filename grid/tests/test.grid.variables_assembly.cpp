#include "boost/hana.hpp"
#include "fub/type_traits.hpp"
#include <type_traits>
#include <typeinfo>

namespace hana = boost::hana;

namespace fub {

template <typename... Ts> using tuple_t = hana::tuple<hana::type<Ts>...>;

template <typename... Ts>
static constexpr tuple_t<Ts...> to_types(const hana::tuple<Ts...>&) noexcept {
  return {};
}

///////////////////////////////////////////////////////////////////////////////
//                                                            [variable_traits]
// {{{

template <typename Collection> struct collection_traits {
private:
  /// This is a helper typedef to detect if `v.as_list()` is a valid expression
  /// for a object v of type V.
  template <typename C>
  using variables_member_function_t = decltype(std::declval<C>().variables());

  // Returns `Variable::as_list()` as it is detected to exist.
  static constexpr auto variables_dispatch(std::true_type,
                                           Collection collection) noexcept {
    return collection.variables();
  }

  /// Returns `hana::tuple<Collection>` because `Collection::variables` does not
  /// exist.
  static constexpr hana::tuple<Collection>
  variables_dispatch(std::false_type, Collection collection) noexcept {
    return {collection};
  }

public:
  /// Returns a `hana::tuple<Vs...>` containing all basic variables which are
  /// contained by this collection.
  static constexpr auto variables(Collection collection) noexcept {
    return variables_dispatch(
        is_detected<variables_member_function_t, Collection>{}, collection);
  }

private:
  /// This is a helper typedef to detect if `v.as_list()` is a valid expression
  /// for a object v of type V.
  template <typename C>
  using size_member_function_t = decltype(std::declval<C>().size());

  // Returns `Variable::as_list()` as it is detected to exist.
  static constexpr auto size_dispatch(std::true_type,
                                      Collection collection) noexcept {
    return collection.size();
  }

  /// Returns `hana::tuple<Collection>` because `Collection::size` does not
  /// exist.
  static constexpr std::ptrdiff_t
  size_dispatch(std::false_type, Collection collection) noexcept {
    return hana::size(variables(collection))();
  }

public:
  /// Returns the number of variables contained by this variable list.
  static constexpr std::ptrdiff_t size(Collection collection) noexcept {
    return size_dispatch(is_detected<size_member_function_t, Collection>{},
                         collection);
  }
};

template <typename> struct is_collection : bool_constant<false> {};

/// Example:
// clang-format off
///     // Construct a simple variable
///     struct Density {};
///     REQUIRE(variable_traits<Density>::as_list() == meta::list<Density>());
///     REQUIRE(variable_traits<Density>::size() == 1);
///
///     // Construct a list containing one variable.
///     struct MyVariableList : basic_variable_list<Density> {};
///     REQUIRE(variable_traits<MyVariableList>::as_list() == meta::list<Density>()); 
///     REQUIRE(variable_traits<MyVariableList>::size() == 1);
// clang-format on
template <typename Variable> class variable_traits {
  /////////////////////////////////////////////////////////////////////////////
  //                                                     [variable_traits.name]
private:
  template <typename V>
  using name_member_function_t = decltype(std::declval<V>().name());

  // Returns Variable::name() because it exists.
  static const char* name_dispatch(std::true_type, Variable variable) noexcept {
    return Variable::name(variable);
  }

  // Returns typeid(Variable).name() because member function `Variable::name()`
  // does NOT exsist.
  static const char* name_dispatch(std::false_type, Variable) noexcept {
    return typeid(Variable).name();
  }

public:
  /// Returns a representative string literal for this variable.
  static const char* name(Variable variable) noexcept {
    return name_dispatch(
        is_detected_exact<const char*, name_member_function_t, Variable>{},
        variable);
  }

private:
  template <typename C, typename V>
  using is_iterable_in_t = decltype(
      std::declval<C>().begin() <= std::declval<C>().next(std::declval<V>()) &&
      std::declval<C>().next(std::declval<V>()) < std::declval<C>().end());

  template <typename Collection>
  static constexpr bool is_iterable_in_dispatch(std::false_type, Collection,
                                                Variable) noexcept {
    return false;
  }

  template <typename Collection>
  static constexpr bool is_iterable_in_dispatch(std::true_type,
                                                Collection collection,
                                                Variable variable) noexcept {
    return collection.begin() <= variable && variable < collection.end();
  }

public:
  template <typename Collection>
  static constexpr bool is_iterable_in(Collection collection,
                                       Variable variable) noexcept {
    return is_iterable_in_dispatch(
        is_detected<is_iterable_in_t, Collection, Variable>{}, collection,
        variable);
  }
};

template <typename... Vs>
struct static_collection : private hana::tuple<Vs...> {
  using Base = hana::tuple<Vs...>;
  using Base::Base;

  static_assert(conjunction<bool_constant<collection_traits<Vs>::size(Vs{}) <=
                                          1>...>::value,
                "Do not merge non-trivial variables.");

  static_assert(!hana::any_of(hana::tuple_t<Vs...>, hana::trait<std::is_const>),
                "Collections should contain non-const variable types.");

  constexpr auto variables() const noexcept {
    const hana::tuple<Vs...>& tuple = *this;
    return tuple;
  }

  /// Note: We need this because there is no constexpr lambd a
  struct plus_size {
    template <typename T>
    constexpr std::size_t operator()(std::size_t n, T collection) const
        noexcept {
      return n + collection_traits<T>::size(collection);
    }
  };

  constexpr std::size_t size() const noexcept {
    return hana::fold(variables(), std::size_t{0}, plus_size{});
  }
};

template <typename... Vs>
struct is_collection<static_collection<Vs...>> : bool_constant<true> {};

template <typename... Vs>
constexpr bool operator==(const static_collection<Vs...>& lhs,
                          const static_collection<Vs...>& rhs) noexcept {
  return lhs.variables() == rhs.variables();
}

template <typename... Vs>
constexpr bool operator!=(const static_collection<Vs...>& lhs,
                          const static_collection<Vs...>& rhs) noexcept {
  return !(lhs == rhs);
}

class dynamic_collection;

class dynamic_iterable_variable {
public:
  dynamic_iterable_variable() = default;

  explicit constexpr dynamic_iterable_variable(std::ptrdiff_t value) noexcept
      : m_value{value} {}

  constexpr std::ptrdiff_t value() const noexcept { return m_value; }

private:
  std::ptrdiff_t m_value;
};

#define FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(op)                  \
  constexpr bool operator op(const dynamic_iterable_variable& lhs,             \
                             const dynamic_iterable_variable& rhs) noexcept {  \
    return lhs.value() op rhs.value();                                         \
  }                                                                            \
  /* */

FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(<);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(<=);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(>);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(>=);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(==);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(!=);

#undef FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE

class dynamic_collection {
public:
  dynamic_collection() = default;

  explicit constexpr dynamic_collection(std::ptrdiff_t size) noexcept
      : m_size{size} {}

  static constexpr auto variables() noexcept {
    return hana::make_tuple(dynamic_iterable_variable{});
  }

  constexpr std::ptrdiff_t size() const noexcept { return m_size; }

  constexpr dynamic_iterable_variable next(dynamic_iterable_variable var) const
      noexcept {
    return dynamic_iterable_variable{var.value() + 1};
  }

  constexpr dynamic_iterable_variable begin() const noexcept {
    return dynamic_iterable_variable{0};
  }

  constexpr dynamic_iterable_variable end() const noexcept {
    return dynamic_iterable_variable{m_size};
  }

public:
  std::ptrdiff_t m_size;
};

constexpr bool operator==(const dynamic_collection& lhs,
                          const dynamic_collection& rhs) noexcept {
  return lhs.size() == rhs.size();
}

constexpr bool operator!=(const dynamic_collection& lhs,
                          const dynamic_collection& rhs) noexcept {
  return lhs.size() != rhs.size();
}

static_assert(
    std::is_nothrow_copy_constructible<dynamic_iterable_variable>::value,
    "dynamic_iterable_variable is not nothrow copy constructible!");

static_assert(is_nothrow_equality_comparable<dynamic_iterable_variable>::value,
              "dynamic_iterable_variable is not equality comparable!");

static_assert(std::is_nothrow_copy_constructible<dynamic_collection>::value,
              "dynamic_collection is not nothrow copy constructible!");
static_assert(is_nothrow_equality_comparable<dynamic_collection>::value,
              "dynamic_collection is not equality comparable!");

static_assert(
    variable_traits<dynamic_iterable_variable>::is_iterable_in(
        dynamic_collection(1), dynamic_iterable_variable(0)),
    "dynamic_iterable_variable is not itrable in dynamic_collection!");

static_assert(
    collection_traits<dynamic_collection>::size(dynamic_collection{}) == 0,
    "collection size does not match in the trait function.");

static_assert(
    collection_traits<dynamic_collection>::size(dynamic_collection(2)) == 2,
    "collection size does not match in the trait function.");

template <typename T> struct collection_comparator {
  T m_value;

  template <typename S>
  constexpr bool equal_to(std::true_type,
                          const collection_comparator<S>& rhs) const noexcept {
    return m_value == rhs.m_value;
  }

  template <typename S>
  constexpr bool equal_to(std::false_type,
                          const collection_comparator<S>& rhs) const noexcept {
    return false;
  }

  template <typename S>
  constexpr bool operator==(const collection_comparator<S>& rhs) const
      noexcept {
    return equal_to(hana::type_c<T> == hana::type_c<S>, rhs);
  }

  template <typename S>
  constexpr bool operator!=(const collection_comparator<S>& rhs) const
      noexcept {
    return !(*this == rhs);
  }
};

struct make_comparator_t {
  template <typename T>
  constexpr collection_comparator<T> operator()(T value) const noexcept {
    return {value};
  }
};
static constexpr make_comparator_t make_comparator{};

/// Returns all available unions as a hana::tuple object.
///
/// Variables are assumed to be light-weight handles to some data. Thus we
/// require some exception guarantees:
///
/// Require: is_nothrow_copy_constructible<Variable>
/// Require: is_nothrow_equality_comparable<const Variable&>
/// @{
template <typename Variable, bool IsCollection>
struct variable_collections_impl : hana::tuple<Variable> {
  using Base = hana::tuple<Variable>;
  using Base::Base;

  static_assert(std::is_nothrow_copy_constructible<Base>{},
                "hana::tuple<Variable> must be nothrow copy constructible.");

  static_assert(is_equality_comparable<Base>{},
                "hana::tuple<Variable> must be nothrow equality comparable.");

  constexpr hana::tuple<Variable> get_collections() const noexcept {
    return *this;
  }
};

template <typename Pred> struct add_one_until_t {
  Pred pred;
  using result_type = std::pair<std::ptrdiff_t, bool>;

  template <typename T>
  constexpr result_type operator()(result_type result, T x) const noexcept {
    std::ptrdiff_t n = result.first;
    bool found = result.second || pred(x);
    if (!found) {
      n += 1;
    }
    return result_type(n, found);
  }
};

template <typename Pred>
constexpr add_one_until_t<Pred> add_one_until(Pred pred) noexcept {
  return {pred};
}

template <typename Xs, typename Pred>
constexpr std::ptrdiff_t index_if(Xs xs, Pred pred) {
  auto x = hana::front(xs);
  if (pred(x)) {
    return 0;
  }
  auto rest = hana::drop_front(xs, int_c<1>);
  return hana::fold_left(rest, std::make_pair(std::ptrdiff_t{1}, false),
                         add_one_until(pred))
      .first;
}

template <template <typename...> class List, typename... Vs>
struct variable_collections_impl<List<Vs...>, false> : hana::tuple<Vs...> {
  using Base = hana::tuple<Vs...>;
  using Base::Base;

  static_assert(std::is_nothrow_copy_constructible<Base>{},
                "hana::tuple<Vs...> must be nothrow copy constructible.");

  static_assert(is_equality_comparable<Base>{},
                "hana::tuple<Vs...> must be nothrow equality comparable.");

  constexpr hana::tuple<Vs...> get_collections() const noexcept {
    return *this;
  }
};

template <typename Variable>
struct variable_collections : variable_collections_impl<Variable, false> {
  using Base = variable_collections_impl<Variable, false>;
  using Base::Base;
};

template <template <typename...> class List, typename... Vs>
struct variable_collections<List<Vs...>>
    : variable_collections_impl<List<Vs...>,
                                is_collection<List<Vs...>>::value> {
  using Base =
      variable_collections_impl<List<Vs...>, is_collection<List<Vs...>>::value>;
  using Base::Base;
};
/// @}

////////////////////////////////////////////////////////////////////////////////
//                                                                 variable_find

/// Given a Variable which is going to be accessed, this class implements how to
/// find the union where the value is stored at. It also computes the index of
/// the union in the set of all present unions.
///
/// Note: This class template works for both cases: For tuple-like types
/// Variables and simple unions.
template <typename Collections>
struct variable_find : variable_collections<Collections> {
  using base = variable_collections<Collections>;
  using base::base;

  // Define Type based Observers

  /// Returns a tuple of types containing all collection types accessible by
  /// this accessor.
  ///
  /// Throws: Nothing.
  static constexpr auto get_collection_types() noexcept {
    return decltype(to_types(
        std::declval<variable_collections<Collections>>().get_collections())){};
  }

  /// Returns an index to the specified collection `C` if it exists in
  /// `Collections`. If it does not exist the function returns `hana::nothing`.
  ///
  /// Throws: Nothing.
  template <typename C>
  static constexpr auto index_collection_type(hana::basic_type<C> c) noexcept {
    return hana::index_if(get_collection_types(), hana::equal.to(c));
  }

  /// This is a unary predicate which checks if a variable is accessible in a
  /// collection.
  template <typename Variable> struct is_accessible_type {
    template <typename Collection>
    constexpr auto operator()(hana::basic_type<Collection> collection) {
      hana::type<Variable> variable;
      return hana::or_(
          collection == variable,
          hana::find(
              to_types(collection_traits<Collection>::variables(Collection{})),
              variable) != hana::nothing);
    }
  };

  /// Returns an optional index if the specified variable type is accessible for
  /// the collection type referred by the returned index.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  static constexpr auto
  index_accessible_collection_type(hana::basic_type<Variable> variable =
                                       hana::basic_type<Variable>()) noexcept {
    return hana::index_if(variable_find::get_collection_types(),
                          variable_find::is_accessible_type<Variable>{});
  }

  template <typename Collection>
  constexpr auto find_collection(Collection collection) const noexcept {
    return hana::find(base::get_collections(), collection);
  }

  /// This is a unary predicate type which checks if a variable is accessible in
  /// a collection.
  template <typename Variable> struct is_accessible {
    collection_comparator<Variable> variable;

    template <typename Collection>
    constexpr bool operator()(Collection c) const noexcept {
      collection_comparator<Collection> collection{c};
      if (collection == variable || variable_traits<Variable>::is_iterable_in(
                                        collection.m_value, variable.m_value)) {
        return true;
      }
      auto vars = collection_traits<Collection>::variables(c);
      if (hana::size(vars) == hana::size_c<1>) {
        auto variables = make_comparator(vars[int_c<0>]);
        if (collection == variables) {
          return false;
        }
      }
      return index_if(vars, is_accessible<Variable>{variable}) <
             collection_traits<Collection>::size(c);
    }
  };

  /// Returns an optional index `I` for which `V` is accessible in
  /// `get_collections()[*I]`.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  constexpr auto index_accessible_collection(Variable variable) const noexcept {
    return index_if(base::get_collections(), is_accessible<Variable>{variable});
  }

  /// Returns an optional union varibale `u` for which `variable` is accessible.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  constexpr auto find_accessible_collection(Variable variable) const noexcept {
    return hana::find_if(base::get_collections(),
                         is_accessible<Variable>{variable});
  }

  template <typename V,
            typename std::enable_if_t<index_accessible_collection_type<V>() !=
                                      hana::nothing>* = nullptr>
  constexpr auto local_index(hana::basic_type<V> var) const noexcept {
    constexpr auto remove_const = hana::metafunction<std::remove_const>;
    auto index = index_accessible_collection(var);
    using Collection = decltype(base::get_collections()[*index]);
    return hana::index_if(collection_traits<Collection>::variables(),
                          hana::equal.to(remove_const(var)));
  }

  /// Returns a tuple of all accessible variables.
  ///
  /// Throws: Nothing.
  constexpr auto get_variables() const noexcept {
    return hana::flatten(
        hana::transform(base::get_collections(), [](auto collection) {
          using C = remove_cvref_t<decltype(collection)>;
          return collection_traits<C>::variables(collection);
        }));
  }

  /// Returns a tuple of all accessible variable types.
  ///
  /// Throws: Nothing.
  constexpr auto get_variable_types() const noexcept {
    return hana::flatten(
        hana::transform(get_collection_types(), [](auto collection) {
          using C = remove_cvref_t<typename decltype(collection)::type>;
          return decltype(
              to_types(collection_traits<C>::variables(std::declval<C>()))){};
        }));
  }
};

struct Density {};
constexpr bool operator==(const Density&, const Density&) noexcept {
  return true;
}
constexpr bool operator!=(const Density&, const Density&) noexcept {
  return false;
}

static_assert(std::is_nothrow_copy_constructible<Density>::value, "");
static_assert(is_nothrow_equality_comparable<Density>::value, "");

struct Momentum {};
constexpr bool operator==(const Momentum&, const Momentum&) noexcept {
  return true;
}
constexpr bool operator!=(const Momentum&, const Momentum&) noexcept {
  return false;
}

static_assert(std::is_nothrow_copy_constructible<Momentum>::value, "");
static_assert(is_nothrow_equality_comparable<Momentum>::value, "");

struct Energy {};
constexpr bool operator==(const Energy&, const Energy&) noexcept {
  return true;
}
constexpr bool operator!=(const Energy&, const Energy&) noexcept {
  return false;
}

static_assert(std::is_nothrow_copy_constructible<Energy>::value, "");
static_assert(is_nothrow_equality_comparable<Energy>::value, "");

struct Pressure {};
constexpr bool operator==(const Pressure&, const Pressure&) noexcept {
  return true;
}
constexpr bool operator!=(const Pressure&, const Pressure&) noexcept {
  return false;
}

static_assert(std::is_nothrow_copy_constructible<Pressure>::value, "");
static_assert(is_nothrow_equality_comparable<Pressure>::value, "");

struct MyCollection : static_collection<Density, dynamic_collection> {
  using base = static_collection<Density, dynamic_collection>;
  using base::base;
};

} // namespace fub

int main() {
  using namespace fub;
  constexpr dynamic_collection collection(8);
  constexpr MyCollection mc(Density{}, collection);
  variable_find<hana::tuple<Pressure, MyCollection>> accessor({}, mc);
  std::ptrdiff_t index =
      accessor.index_accessible_collection(dynamic_iterable_variable(8));
  return index;
}
