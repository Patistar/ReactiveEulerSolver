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

#include <vector>

namespace fub {
inline namespace v1 {
namespace p4est {
namespace {
std::ptrdiff_t local_index(const ghost_layer<2>& ghost, int which_tree,
                           int local_num) noexcept {
  span<const int> offsets = ghost.quadrants_tree_offsets();
  assert(which_tree + 1 < offsets.size());
  return offsets[which_tree] + local_num;
}

std::ptrdiff_t local_index(const p4est_ghost_t& ghost, int which_tree,
                           int local_num) noexcept {
  std::ptrdiff_t offset = ghost.tree_offsets[which_tree];
  return offset + local_num;
}

void add_neighbors(span<mesh<2>::face_neighbors_t> quad_to_face_neighbors,
                   const p4est_ghost_t& ghost_layer, bool lhs_is_ghost,
                   const quadrant<2>& lhs, int lhs_face, bool rhs_is_ghost,
                   const quadrant<2>& rhs, int rhs_face) {
  if (!lhs_is_ghost) {
    std::ptrdiff_t lhs_idx =
        local_index(ghost_layer, lhs.which_tree(), lhs.local_num());
    quad_to_face_neighbors[lhs_idx][lhs_face].push_back(rhs);
  }
  if (!rhs_is_ghost) {
    std::ptrdiff_t rhs_idx =
        local_index(ghost_layer, rhs.which_tree(), rhs.local_num());
    quad_to_face_neighbors[rhs_idx][rhs_face].push_back(lhs);
  }
}

void add_neighbors(span<mesh<2>::face_neighbors_t> quad_to_face_neighbors,
                   p4est_iter_face_info_t& info) {
  p4est_iter_face_side_t& left = *p4est_iter_fside_array_index(&info.sides, 0);
  p4est_iter_face_side_t& right = *p4est_iter_fside_array_index(&info.sides, 1);
  if (left.is_hanging) {
    add_neighbors(quad_to_face_neighbors, *info.ghost_layer,
                  left.is.hanging.is_ghost[0], *left.is.hanging.quad[0],
                  left.face, right.is.full.is_ghost, *right.is.full.quad,
                  right.face);
    add_neighbors(quad_to_face_neighbors, *info.ghost_layer,
                  left.is.hanging.is_ghost[1], *left.is.hanging.quad[1],
                  left.face, right.is.full.is_ghost, *right.is.full.quad,
                  right.face);
  } else if (right.is_hanging) {
    add_neighbors(quad_to_face_neighbors, *info.ghost_layer,
                  right.is.hanging.is_ghost[0], *right.is.hanging.quad[0],
                  right.face, left.is.full.is_ghost, *left.is.full.quad,
                  left.face);
    add_neighbors(quad_to_face_neighbors, *info.ghost_layer,
                  right.is.hanging.is_ghost[1], *right.is.hanging.quad[1],
                  right.face, left.is.full.is_ghost, *left.is.full.quad,
                  left.face);
  } else {
    add_neighbors(quad_to_face_neighbors, *info.ghost_layer,
                  left.is.full.is_ghost, *left.is.full.quad, left.face,
                  right.is.full.is_ghost, *right.is.full.quad, right.face);
  }
}
} // namespace

mesh<2>::mesh(fub::p4est::forest<2> forest, fub::p4est::ghost_layer<2> ghost,
              mesh_data_size size, allocator_type alloc)
    : m_forest{std::move(forest)}, m_ghost{std::move(ghost)}, m_allocator{
                                                                  alloc} {
  allocate(size, alloc);
}

mesh<2>::mesh(mesh&& other) noexcept
    : m_forest{std::move(other.m_forest)}, m_ghost{std::move(other.m_ghost)},
      m_allocator{other.m_allocator}, m_local_to_data{std::exchange(
                                          other.m_local_to_data, {})},
      m_ghost_to_data{std::exchange(other.m_ghost_to_data, {})},
      m_mirror_to_data{std::exchange(other.m_mirror_to_data, {})},
      m_quad_to_face_neighbors{
          std::exchange(other.m_quad_to_face_neighbors, {})},
      m_level_to_local{std::exchange(other.m_level_to_local, {})},
      m_coarse_fine_interfaces{
          std::exchange(other.m_coarse_fine_interfaces, {})} {}

mesh<2>& mesh<2>::operator=(mesh&& other) noexcept {
  deallocate();
  m_forest = std::move(other.m_forest);
  m_ghost = std::move(other.m_ghost);
  m_allocator = other.m_allocator;
  m_local_to_data = std::exchange(other.m_local_to_data, {});
  m_ghost_to_data = std::exchange(other.m_ghost_to_data, {});
  m_mirror_to_data = std::exchange(other.m_mirror_to_data, {});
  m_quad_to_face_neighbors = std::exchange(other.m_quad_to_face_neighbors, {});
  m_level_to_local = std::exchange(other.m_level_to_local, {});
  m_coarse_fine_interfaces = std::exchange(other.m_coarse_fine_interfaces, {});
  return *this;
}

mesh<2>::~mesh() noexcept { deallocate(); }

namespace {
template <typename A, typename T>
using rebind_alloc =
    typename std::allocator_traits<A>::template rebind_alloc<T>;

template <typename T, typename Allocator>
span<T> allocate_(Allocator alloc, std::ptrdiff_t size) {
  rebind_alloc<Allocator, T> a(alloc);
  return span<T>(a.allocate(size), size);
}

template <typename Allocator, typename T>
void deallocate_(Allocator alloc, fub::span<T> span) noexcept {
  if (span) {
    rebind_alloc<Allocator, T> a(alloc);
    a.deallocate(span.data(), span.size());
  }
}
} // namespace

void mesh<2>::deallocate() noexcept {
  deallocate_(m_allocator, m_local_to_data);
  deallocate_(m_allocator, m_ghost_to_data);
  deallocate_(m_allocator, m_mirror_to_data);
  deallocate_(m_allocator, m_quad_to_face_neighbors);
  for (span<std::ptrdiff_t> level : m_level_to_local) {
    deallocate_(m_allocator, level);
  }
  deallocate_(m_allocator, m_coarse_fine_interfaces);
}

void mesh<2>::allocate(mesh_data_size size, allocator_type alloc) {
  span<span<byte>> local_to_data{};
  span<span<byte>> ghost_to_data{};
  span<span<byte>> mirror_to_data{};
  span<face_neighbors_t> quad_to_face_neighbors{};
  std::array<span<std::ptrdiff_t>, P4EST_QMAXLEVEL> level_to_local{};
  span<face_info<2>> coarse_fine_interfaces{};
  try {
    local_to_data =
        allocate_<span<byte>>(alloc, m_forest.local_num_quadrants());
    for (span<byte>& data : local_to_data) {
      data = allocate_<byte>(alloc, size.local);
    }
    ghost_to_data = allocate_<span<byte>>(alloc, m_ghost.quadrants().size());
    for (span<byte>& data : ghost_to_data) {
      data = allocate_<byte>(alloc, size.ghost);
    }
    mirror_to_data = allocate_<span<byte>>(alloc, m_ghost.mirrors().size());
    for (span<byte>& data : mirror_to_data) {
      data = allocate_<byte>(alloc, size.ghost);
    }
    quad_to_face_neighbors =
        allocate_<face_neighbors_t>(alloc, m_forest.local_num_quadrants());
    new (quad_to_face_neighbors.data())
        face_neighbors_t[m_forest.local_num_quadrants()];
    struct counter_context {
      std::array<std::ptrdiff_t, P4EST_QMAXLEVEL> level_counter;
      std::ptrdiff_t hanging_counter;
    };
    counter_context ctx{};
    p4est_iterate(m_forest.native(), m_ghost.native(), &ctx,
                  [](p4est_iter_volume_info_t* info, void* ctx_ptr) {
                    auto ctx = static_cast<counter_context*>(ctx_ptr);
                    info->quad->p.piggy3.which_tree = info->treeid;
                    info->quad->p.piggy3.local_num = info->quadid;
                    quadrant<2> quad{*info->quad};
                    ctx->level_counter[quad.level()] += 1;
                  },
                  [](p4est_iter_face_info_t* info, void* ctx_ptr) {
                    auto ctx = static_cast<counter_context*>(ctx_ptr);
                    if (!info->tree_boundary) {
                      p4est_iter_face_side_t* left =
                          p4est_iter_fside_array_index(&info->sides, 0);
                      p4est_iter_face_side_t* right =
                          p4est_iter_fside_array_index(&info->sides, 1);
                      if (left->is_hanging || right->is_hanging) {
                        ctx->hanging_counter += 1;
                      }
                    }
                  },
                  nullptr);
    for (int level = 0; level < P4EST_QMAXLEVEL; ++level) {
      level_to_local[level] =
          allocate_<std::ptrdiff_t>(alloc, ctx.level_counter[level]);
    }
    coarse_fine_interfaces =
        allocate_<face_info<2>>(alloc, ctx.hanging_counter);
    struct fill_context {
      const fub::p4est::ghost_layer<2>* ghost;
      std::array<span<std::ptrdiff_t>, P4EST_QMAXLEVEL>* level_to_local;
      std::array<std::ptrdiff_t, P4EST_QMAXLEVEL> level_counter;
      span<face_info<2>> coarse_fine_interfaces;
      std::ptrdiff_t interface_counter;
      span<face_neighbors_t> quad_to_face_neighbors;
    };
    fill_context fill_ctx{&m_ghost, &level_to_local,
                          {0},      coarse_fine_interfaces,
                          0,        quad_to_face_neighbors};
    p4est_iterate(
        m_forest.native(), m_ghost.native(), &fill_ctx,
        [](p4est_iter_volume_info_t* info, void* ctx_ptr) {
          auto ctx = static_cast<fill_context*>(ctx_ptr);
          const quadrant<2> quad{*info->quad};
          const int level = quad.level();
          const int i = ctx->level_counter[level];
          (*ctx->level_to_local)[level][i] =
              local_index(*ctx->ghost, quad.which_tree(), quad.local_num());
          ctx->level_counter[level] += 1;
        },
        [](p4est_iter_face_info_t* info, void* ctx_ptr) {
          auto ctx = static_cast<fill_context*>(ctx_ptr);
          if (!info->tree_boundary) {
            add_neighbors(ctx->quad_to_face_neighbors, *info);
            p4est_iter_face_side_t* left =
                p4est_iter_fside_array_index(&info->sides, 0);
            p4est_iter_face_side_t* right =
                p4est_iter_fside_array_index(&info->sides, 1);
            if (right->is_hanging) {
              std::swap(left, right);
            }
            if (left->is_hanging) {
              const std::ptrdiff_t i = ctx->interface_counter;
              // if left is hanging, right can not be.
              assert(!right->is_hanging);
              face_info_data<2> coarse{*right->is.full.quad,
                                       bool(right->is.full.is_ghost)};
              std::array<face_info_data<2>, 2> fine{
                  face_info_data<2>{*left->is.hanging.quad[0],
                                    bool(left->is.hanging.is_ghost[0])},
                  face_info_data<2>{*left->is.hanging.quad[1],
                                    bool(left->is.hanging.is_ghost[1])}};
              face_info<2> info(coarse, fine, fub::face(right->face));
              ctx->coarse_fine_interfaces[i] = info;
              ctx->interface_counter += 1;
              // Add neighbor relationship
            }
          }
        },
        nullptr);
    deallocate();
    m_local_to_data = local_to_data;
    m_ghost_to_data = ghost_to_data;
    m_mirror_to_data = mirror_to_data;
    m_quad_to_face_neighbors = quad_to_face_neighbors;
    m_level_to_local = level_to_local;
    m_coarse_fine_interfaces = coarse_fine_interfaces;
    m_allocator = alloc;
  } catch (...) {
    deallocate_(alloc, local_to_data);
    deallocate_(alloc, ghost_to_data);
    deallocate_(alloc, mirror_to_data);
    deallocate_(alloc, quad_to_face_neighbors);
    for (span<std::ptrdiff_t> level : level_to_local) {
      deallocate_(alloc, level);
    }
    deallocate_(alloc, coarse_fine_interfaces);
  }
}

void mesh<2>::reset(fub::p4est::forest<2> forest,
                    fub::p4est::ghost_layer<2> ghost, mesh_data_size size,
                    allocator_type alloc) {
  m_forest = std::move(forest);
  m_ghost = std::move(ghost);
  allocate(size, alloc);
}

namespace {
void mpi_isend(span<const byte> span, int dest, int tag, MPI_Comm comm,
               std::vector<MPI_Request>& reqs) {
  MPI_Request req{};
  MPI_Isend(span.data(), span.size(), MPI_BYTE, dest, tag, comm, &req);
  reqs.push_back(req);
}

void mpi_irecv(span<byte> span, int dest, int tag, MPI_Comm comm,
               std::vector<MPI_Request>& reqs) {
  MPI_Request req{};
  MPI_Irecv(span.data(), span.size(), MPI_BYTE, dest, tag, comm, &req);
  reqs.push_back(req);
}
} // namespace

int mesh<2>::synchronize_ghost_layer(projection_type projection) {
  auto mirrors = m_ghost.mirrors();
  auto mirroridx_by_process = m_ghost.mirrors_by_process();
  auto mirror_offsets = m_ghost.mirrors_by_process_offsets();
  std::vector<MPI_Request> requests;
  requests.reserve(mirroridx_by_process.size() + m_ghost.quadrants().size());
  for (int rank = 0; rank + 1 < mirror_offsets.size(); ++rank) {
    int lower = mirror_offsets[rank];
    int upper = mirror_offsets[rank + 1];
    for (int i = lower; i < upper; ++i) {
      int mirroridx = mirroridx_by_process[i];
      quadrant<2> mirror = mirrors[mirroridx];
      optional<std::ptrdiff_t> locidx = find(forest(), mirror);
      assert(locidx);
      fub::invoke(projection, mirror, m_local_to_data[*locidx],
                  m_mirror_to_data[mirroridx]);
      mpi_isend(m_mirror_to_data[mirroridx], rank,
                mirrors[mirroridx].local_num(), forest().mpi_communicator(),
                requests);
    }
  }
  span<const quadrant<2>> ghost = m_ghost.quadrants();
  span<const int> ghost_offsets = m_ghost.quadrants_process_offsets();
  for (int rank = 0; rank + 1 < ghost_offsets.size(); ++rank) {
    int lower = ghost_offsets[rank];
    int upper = ghost_offsets[rank + 1];
    for (int i = lower; i < upper; ++i) {
      mpi_irecv(m_ghost_to_data[i], rank, ghost[i].local_num(),
                forest().mpi_communicator(), requests);
    }
  }
  return MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

const forest<2>& mesh<2>::forest() const noexcept { return m_forest; }

const ghost_layer<2>& mesh<2>::ghost_layer() const noexcept { return m_ghost; }

span<byte> mesh<2>::local_data(quadrant<2> quad) noexcept {
  std::ptrdiff_t idx =
      local_index(m_ghost, quad.which_tree(), quad.local_num());
  assert(0 <= idx && idx < m_local_to_data.size());
  return m_local_to_data[idx];
}

span<const byte> mesh<2>::local_data(quadrant<2> quad) const noexcept {
  std::ptrdiff_t idx =
      local_index(m_ghost, quad.which_tree(), quad.local_num());
  assert(0 <= idx && idx < m_local_to_data.size());
  return m_local_to_data[idx];
}

span<byte> mesh<2>::local_data(std::ptrdiff_t idx) noexcept {
  assert(0 <= idx && idx < m_local_to_data.size());
  return m_local_to_data[idx];
}

span<const byte> mesh<2>::local_data(std::ptrdiff_t idx) const noexcept {
  assert(0 <= idx && idx < m_local_to_data.size());
  return m_local_to_data[idx];
}

span<byte> mesh<2>::ghost_data(quadrant<2> quad) noexcept {
  std::ptrdiff_t idx = quad.local_num();
  assert(0 <= idx && idx < m_ghost_to_data.size());
  return m_ghost_to_data[idx];
}

span<const byte> mesh<2>::ghost_data(quadrant<2> quad) const noexcept {
  std::ptrdiff_t idx = quad.local_num();
  assert(0 <= idx && idx < m_ghost_to_data.size());
  return m_ghost_to_data[idx];
}

span<const quadrant<2>> mesh<2>::face_neighbors(quadrant<2> quad, face f) const
    noexcept {
  const std::ptrdiff_t idx =
      local_index(m_ghost, quad.which_tree(), quad.local_num());
  const span<const quadrant<2>> span(m_quad_to_face_neighbors[idx][f]);
  return span;
}

span<const std::ptrdiff_t> mesh<2>::quadrants_at_level(int level) const
    noexcept {
  return m_level_to_local[level];
}

span<const face_info<2>> mesh<2>::coarse_fine_interfaces() const noexcept {
  return m_coarse_fine_interfaces;
}
} // namespace p4est
} // namespace v1
} // namespace fub