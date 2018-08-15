#include "fub/p4est/quadrant.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {
namespace {
uint32_t reverse(uint32_t x, int bits) noexcept {
  x = ((x & 0x55555555) << 1) | ((x & 0xAAAAAAAA) >> 1);
  x = ((x & 0x33333333) << 2) | ((x & 0xCCCCCCCC) >> 2);
  x = ((x & 0x0F0F0F0F) << 4) | ((x & 0xF0F0F0F0) >> 4);
  x = ((x & 0x00FF00FF) << 8) | ((x & 0xFF00FF00) >> 8);
  x = ((x & 0x0000FFFF) << 16) | ((x & 0xFFFF0000) >> 16);
  return x >> (32 - bits);
}
} // namespace

int quadrant<2>::which_tree() const noexcept { return m_native.p.which_tree; }

int quadrant<2>::local_num() const noexcept {
  return m_native.p.piggy3.local_num;
}

int quadrant<2>::level() const noexcept { return m_native.level; }

std::array<int, 2> quadrant<2>::coordinates() const noexcept {
  int x = reverse(reverse(m_native.x, 30), level());
  int y = reverse(reverse(m_native.y, 30), level());
  return {x, y};
}

bool operator==(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept {
  return p4est_quadrant_is_equal(&lhs.native(), &rhs.native());
}

bool operator!=(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept {
  return !(lhs == rhs);
}

bool operator<(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept {
  return p4est_quadrant_compare(&lhs.native(), &rhs.native()) < 0;
}

quadrant<2> parent(const quadrant<2>& quad) noexcept {
  p4est_quadrant_t p;
  p4est_quadrant_parent(&quad.native(), &p);
  return p;
}

quadrant<2> face_neighbor(const quadrant<2>& quad, int face) noexcept {
  p4est_quadrant_t nb;
  p4est_quadrant_face_neighbor(&quad.native(), face, &nb);
  return quadrant<2>(nb);
}

optional<face> find_adjacent_face(const quadrant<2>& left,
                                  const quadrant<2>& right) noexcept {
  p4est_quadrant_t l = left.native();
  p4est_quadrant_t r = right.native();
  if (right.level() < left.level()) {
    p4est_quadrant_ancestor(&left.native(), right.level(), &l);
  } else if (left.level() < right.level()) {
    p4est_quadrant_ancestor(&right.native(), left.level(), &r);
  }
  quadrant<2> lhs(l);
  quadrant<2> rhs(r);
  for (int dim = 0; dim < 2; ++dim) {
    for (int dir = 0; dir < 2; ++dir) {
      face f = face{axis(dim), direction(dir)};
      if (face_neighbor(lhs, f) == rhs) {
        return f;
      }
    }
  }
  return nullopt;
}

int child_id(const quadrant<2>& quad) noexcept {
  return p4est_quadrant_child_id(&quad.native());
}

std::array<quadrant<2>, 4> children(const quadrant<2>& quad) noexcept {
  std::array<p4est_quadrant_t, 4> natives{};
  p4est_quadrant_childrenv(&quad.native(), natives.data());
  for (p4est_quadrant_t& child : natives) {
    child.p.piggy3.which_tree = quad.which_tree();
  }
  std::array<quadrant<2>, 4> array{natives[0], natives[1], natives[2],
                                   natives[3]};
  return array;
}

} // namespace p4est
} // namespace v1
} // namespace fub