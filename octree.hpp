#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "morton.hpp"
#include <map>
#include <array>
#include <vector>
#include <iostream>

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
    using iterator = typename tree_type::iterator;
    using const_iterator = typename tree_type::const_iterator;

    Octree() { initializeMasks(); }
    const array<Morton<Dim, N>, N> &getMasksRightTrue() const { return masks_right_true_; }
    const array<Morton<Dim, N>, N> &getMasksRightFalse() const { return masks_right_false_; }
    const array<Morton<Dim, N>, N> &getMasksSelectLevel() const { return masks_select_level_; }
    const array<Morton<Dim, N>, N> &getMasksUnselectLevel() const { return masks_unselect_level_; }
    const array<Morton<Dim, N>, Dim> &getMasksDim() const { return masks_dim_; }
    const array<array<Morton<Dim, N>, numOctants>,N> &getOctants() const { return octants_;}
    bool insert(const value_type &value, size_t level);
    size_t getLevel(const Morton<Dim, N> &m) const{
        //logic for a complete tree a cell is always part of an octant
        size_t level = m.countr_zero();
        if (level==0)
            return level;
        size_t subLevel = level-1;
        Morton<Dim, N> mSub = 1;
        mSub <<= subLevel;
        while (subLevel != -1)
        {
            if (tree_.count(m+mSub) ==0){
                return level;
            }
            level = subLevel;
            subLevel -= 1;
            mSub >>= 1;
        }
        return level;
    }

    void print(){
        for(auto itr = tree_.cbegin(); itr != tree_.cend();++itr)
            std::cout << itr->first << " " << itr->second << std::endl;
    }

    auto getCells() const{
        std::vector<std::array< std::array<int,Dim> ,2>> cells;
        cells.reserve(tree_.size());
        std::array<std::array<int,Dim>,2> cell;
        for(auto itr = tree_.cbegin(); itr != tree_.cend();++itr){
            auto level = getLevel(itr->first);
            cell[0][0] = encoder_.template decode<0>(itr->first).template getIntValue<int>();
            cell[0][1] = encoder_.template decode<1>(itr->first).template getIntValue<int>();
            cell[1][0] = cell[0][0] + (1<<level);
            cell[1][1] = cell[0][1] + (1<<level);
            cells.push_back(cell);
        }
    return cells;
    }
private:
    array<Morton<Dim, N>, N> masks_right_true_;
    array<Morton<Dim, N>, N> masks_right_false_;
    array<Morton<Dim, N>, N> masks_select_level_;
    array<Morton<Dim, N>, N> masks_unselect_level_;
    array<Morton<Dim, N>, Dim> masks_dim_;
    array<array<Morton<Dim, N>, numOctants>, N> octants_;
    MortonEncoder<Dim, N> encoder_;
    void initializeMasks();
    tree_type tree_;
};

template <size_t Dim, size_t N, typename T>
void Octree<Dim, N, T>::initializeMasks()
{
    auto ones = static_cast<Morton<Dim,N> >((static_cast<size_t>(1) << Dim)-1);
    Morton<Dim,N> morton;
    for (auto i(0); i < N; ++i)
    {
        morton |= ones;
        masks_right_true_[i] = morton;
        masks_right_false_[i] = ~morton;
        morton <<= 1;
    }
    morton = ones;
    for (auto i(0); i < N; ++i)
    {
        masks_select_level_[i] = morton;
        masks_unselect_level_[i] = ~morton;
        morton <<= 1;
    }
    morton = 0; morton -= 1;
    for (auto i(0); i < Dim; ++i)
    {
        masks_dim_[i] = morton;
        morton.BitArrayBase<Morton<Dim, N>, Dim * N>::operator<<=(1); 
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
    Morton<Dim, N> mortonNew(value.first), mortonLoc(value.first);
    mortonNew &=  masks_right_false_[level];
    mortonLoc.set(0, level); // put a one to be sure of the level
    size_t levelUp = N-1;
    auto itr = tree_.begin();
    if (!tree_.empty())
    {
        itr = tree_.upper_bound(mortonLoc);
        Morton<Dim, N> mortonParent=(--itr)->first;
        for(levelUp = level; levelUp<N && (mortonLoc != mortonParent); ++levelUp){
            mortonLoc &= masks_right_false_[levelUp];
        }
    } //levelUp is the octant of the highest level that contains the to be inserted cell
    mortonLoc = mortonNew;
    for (; levelUp != (level-1) ; --levelUp){
        for (auto i(0); i< numOctants; ++i){
            mortonLoc &= masks_right_false_[levelUp];
            mortonLoc |= octants_[levelUp][i];
            itr = tree_.insert(itr, {mortonLoc, value.second});
        }
    }
    return true;
};

#endif // OCTREE_HPP