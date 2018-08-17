// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/p4est/mesh.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {
namespace {
std::ptrdiff_t local_index(const p4est_t& forest, int treeid, int quadid) {
  return p4est_tree_array_index(forest.trees, treeid)->quadrants_offset +
         quadid;
}

void add_neighbors(span<mesh<2>::face_neighbors_t> quad_to_face_neighbors,
                   const p4est_t& forest, bool lhs_is_ghost,
                   const quadrant<2>& lhs, int lhs_face, bool rhs_is_ghost,
                   const quadrant<2>& rhs, int rhs_face) {
  if (!rhs_is_ghost || !lhs_is_ghost) {
    if (rhs_is_ghost) {
      std::ptrdiff_t lhs_idx =
          local_index(forest, lhs.which_tree(), lhs.local_num());
      quad_to_face_neighbors[lhs_idx][lhs_face].push_back({rhs, rhs_is_ghost});
    } else {
      std::ptrdiff_t rhs_idx =
          local_index(forest, rhs.which_tree(), rhs.local_num());
      quad_to_face_neighbors[rhs_idx][rhs_face].push_back({lhs, lhs_is_ghost});
    }
  }
}

void add_neighbors(span<mesh<2>::face_neighbors_t> quad_to_face_neighbors,
                   p4est_iter_face_info_t& info) {
  p4est_iter_face_side_t& left = *p4est_iter_fside_array_index(&info.sides, 0);
  p4est_iter_face_side_t& right = *p4est_iter_fside_array_index(&info.sides, 1);
  if (left.is_hanging) {
    add_neighbors(quad_to_face_neighbors, *info.p4est,
                  left.is.hanging.is_ghost[0], *left.is.hanging.quad[0],
                  left.face, right.is.full.is_ghost, *right.is.full.quad,
                  right.face);
    add_neighbors(quad_to_face_neighbors, *info.p4est,
                  left.is.hanging.is_ghost[1], *left.is.hanging.quad[1],
                  left.face, right.is.full.is_ghost, *right.is.full.quad,
                  right.face);
  } else if (right.is_hanging) {
    add_neighbors(quad_to_face_neighbors, *info.p4est,
                  right.is.hanging.is_ghost[0], *right.is.hanging.quad[0],
                  right.face, left.is.full.is_ghost, *left.is.full.quad,
                  left.face);
    add_neighbors(quad_to_face_neighbors, *info.p4est,
                  right.is.hanging.is_ghost[1], *right.is.hanging.quad[1],
                  right.face, left.is.full.is_ghost, *left.is.full.quad,
                  left.face);
  } else {
    add_neighbors(quad_to_face_neighbors, *info.p4est, left.is.full.is_ghost,
                  *left.is.full.quad, left.face, right.is.full.is_ghost,
                  *right.is.full.quad, right.face);
  }
}

std::vector<mesh<2>::face_neighbors_t>
make_face_neighbors(fub::v1::p4est::forest<2>& forest,
                    fub::v1::p4est::ghost_layer<2>& ghost_layer) {
  std::vector<mesh<2>::face_neighbors_t> quad_to_face_neighbors(
      forest.local_num_quadrants());
  p4est_iterate(forest.native_handle(), ghost_layer.native_handle(),
                &quad_to_face_neighbors, nullptr,
                [](p4est_iter_face_info_t* info, void* ctx_ptr) {
                  auto quad_to_face_neighbors =
                      static_cast<std::vector<mesh<2>::face_neighbors_t>*>(
                          ctx_ptr);
                  if (!info->tree_boundary) {
                    add_neighbors(*quad_to_face_neighbors, *info);
                  }
                },
                nullptr);
  return quad_to_face_neighbors;
}

std::vector<coarse_fine_interface<2>>
make_coarse_fine_interfaces(fub::v1::p4est::forest<2>& forest,
                            fub::v1::p4est::ghost_layer<2>& ghost_layer) {
  std::vector<coarse_fine_interface<2>> coarse_fine_interfaces;
  p4est_iterate(
      forest.native_handle(), ghost_layer.native_handle(),
      &coarse_fine_interfaces, nullptr,
      [](p4est_iter_face_info_t* info, void* ctx_ptr) {
        auto coarse_fine_interfaces =
            static_cast<std::vector<coarse_fine_interface<2>>*>(ctx_ptr);
        if (!info->tree_boundary) {
          p4est_iter_face_side_t* left =
              p4est_iter_fside_array_index(&info->sides, 0);
          p4est_iter_face_side_t* right =
              p4est_iter_fside_array_index(&info->sides, 1);
          if (right->is_hanging) {
            std::swap(left, right);
          }
          if (left->is_hanging) {
            // if left is hanging, right can not be.
            assert(!right->is_hanging);
            mesh_quadrant<2> coarse{*right->is.full.quad,
                                    bool(right->is.full.is_ghost)};
            std::array<mesh_quadrant<2>, 2> fine{
                mesh_quadrant<2>{*left->is.hanging.quad[0],
                                 bool(left->is.hanging.is_ghost[0])},
                mesh_quadrant<2>{*left->is.hanging.quad[1],
                                 bool(left->is.hanging.is_ghost[1])}};
            coarse_fine_interfaces->push_back(
                coarse_fine_interface<2>{{fub::face(right->face), coarse},
                                         {fub::face(left->face), fine}});
          }
        }
      },
      nullptr);
  coarse_fine_interfaces.shrink_to_fit();
  return coarse_fine_interfaces;
}
} // namespace

mesh<2>::mesh(fub::p4est::forest<2>& forest, fub::p4est::ghost_layer<2>& ghost)
    : m_quad_to_face_neighbors{make_face_neighbors(forest, ghost)},
      m_coarse_fine_interfaces{make_coarse_fine_interfaces(forest, ghost)} {}

} // namespace p4est
} // namespace v1
} // namespace fub