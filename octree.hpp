#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "morton.hpp"
#include <map>
#include <array>

using std::array;
using std::map;

template <size_t Dim, size_t N, typename T>
class Octree
{
public:
    constexpr static size_t numOctants = (static_cast<size_t>(1) << Dim);
    using tree_type = map<Morton<Dim, N>, T>;
    using key_type = typename tree_type::key_type;
    using mapped_type = typename tree_type::mapped_type;
    using value_type = typename tree_type::value_type;
    Octree() { initializeMasks(); }
    const array<Morton<Dim, N>, N> &getMasksRightTrue() const { return masks_right_true_; }
    const array<Morton<Dim, N>, N> &getMasksLeftTrue() const { return masks_left_true_; }
    const array<Morton<Dim, N>, N> &getMasksSelectLevel() const { return masks_select_level_; }
    const array<Morton<Dim, N>, N> &getMasksUnselectLevel() const { return masks_unselect_level_; }
    const array<Morton<Dim, N>, Dim> &getMasksDim() const { return masks_dim_; }
    const array<array<Morton<Dim, N>, numOctants>,N> &getOctants() const { return octants_;}
    bool insert(const value_type &value, size_t level);
private:
    array<Morton<Dim, N>, N> masks_right_true_;
    array<Morton<Dim, N>, N> masks_left_true_;
    array<Morton<Dim, N>, N> masks_select_level_;
    array<Morton<Dim, N>, N> masks_unselect_level_;
    array<Morton<Dim, N>, Dim> masks_dim_;
    array<array<Morton<Dim, N>, numOctants>, N> octants_;

    void initializeMasks();
    tree_type tree_;
};

template <size_t Dim, size_t N, typename T>
void Octree<Dim, N, T>::initializeMasks()
{
    Morton<Dim, N> morton(Morton<Dim, N>::ones);
    for (auto i(0); i < N; ++i)
    {
        morton |= Morton<Dim, N>::ones;
        masks_right_true_[i] = morton;
        masks_left_true_[i] = ~morton;
        morton <<= 1;
    }
    morton = Morton<Dim, N>::ones;
    for (auto i(0); i < N; ++i)
    {
        masks_select_level_[i] = morton;
        masks_unselect_level_[i] = ~morton;
        morton <<= 1;
    }
    morton = 0;
    morton -= 1;
    for (auto i(0); i < Dim; ++i)
    {
        masks_dim_[i] = morton;
        morton.bitset<Dim*N>::operator<<=(1); 
    }
    array<Morton<Dim, N>, numOctants> octants;
    for (auto i(0); i < numOctants; ++i)
        octants[i] = Morton<Dim, N>(i);
    for (auto j(0); j < N; ++j)
    {
        octants_[j] = octants;
        for (auto i(0); i < numOctants; ++i)
            octants[i] <<= 1;
    }
};

template <size_t Dim, size_t N, typename T>
bool Octree<Dim, N, T>::insert(const value_type &value, size_t level)
{
    Morton<Dim, N> morton;
    morton = value.first;
    morton.set(0, level); // put a one to be sure of the level
    if (tree_.contains(morton))
        return false;
    for (auto i(level); i < N; ++i)
    {
        morton &= masks_left_true_[level];
        for (auto j(0); j < numOctants; ++j)
        {
            Morton<Dim, N> morton_oct = (morton | octants_[i][j]);
            tree_.insert({morton_oct, value.second});
        }
        morton &= masks_left_true_[level];
    }

    return true;
};

#endif // OCTREE_HPP