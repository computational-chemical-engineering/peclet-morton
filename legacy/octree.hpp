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
    struct Cell
    {
        int level;
        T* ptr;
    };
    constexpr static size_t numOctants = (static_cast<size_t>(1) << Dim);
    using cell_type = typename Octree::Cell;
    using tree_type = map<Morton<Dim, N>, Cell>;
    using key_type = typename tree_type::key_type;
    using mapped_type = typename tree_type::mapped_type;
    using value_type = typename tree_type::value_type;
    using iterator = typename tree_type::iterator;
    using const_iterator = typename tree_type::const_iterator;

    Octree():levelCoarse_(N) { initializeMasks(); }
//    const array<Morton<Dim, N>, N> &getMasksRightTrue() const { return masks_right_true_; }
    const array<Morton<Dim, N>, N> &getMasksGreaterThan() const { return masks_greater_than_; }
    const array<Morton<Dim, N>, N> &getMasksGreaterEqual() const { return masks_greater_equal_; }
    const array<Morton<Dim, N>, N> &getMasksSelectLevel() const { return masks_select_level_; }
    const array<Morton<Dim, N>, N> &getMasksUnselectLevel() const { return masks_unselect_level_; }
    const array<array<Morton<Dim, N>, N>,Dim> &getMasksUnits() const { return masks_units_; }
    const array<Morton<Dim, N>, Dim> &getMasksDim() const { return masks_dim_; }
    const array<array<Morton<Dim, N>, numOctants>,N> &getOctants() const { return octants_;}
    const MortonEncoder<Dim, N> &getMortonEncoder() const {return encoder_;}

    void setCoarseGrid(size_t level, array<size_t, Dim> numCells)
    {
        levelCoarse_ = level;
        size_t maxNumCells = (1 << (N-level));
        array<Morton<Dim,N>, Dim> dx;
        int unit = 1;
        for(auto i(0); i<Dim; ++i)
        {
            if (numCells[i] >= maxNumCells)
            {
                throw std::out_of_range("Cannot accomodate the required coarse cells. Increase the number of bits for the morton numbers.");
            }
            dx[i] = unit;
            dx[i] <<= level;
            unit <<= 1;
        }
        Cell cell;
        cell.level = level;
        cell.ptr = NULL;
        Morton<Dim,N> m;
        size_t numCells2(1);
        Morton<Dim, N> dx2;
        if constexpr (Dim==3)
        {
            numCells2 = numCells[2];
            auto dx2 = dx[2];
        }
        auto itr = tree_.begin();
        for(auto i2=0; i2<numCells2; ++i2, m += dx2){
            m &= ~masks_dim_[1];
            for(auto i1=0; i1<numCells[1]; ++i1, m += dx[1]){
                m &= ~masks_dim_[0];
                for(auto i0=0; i0<numCells[0]; ++i0, m += dx[0]){
                    itr = tree_.insert(itr, {m,cell});
                    //std::cout << "inserted: " << m << " " << cell.level << std::endl;

                }
            }
        }
    }

    void findCellsInParent(const Morton<Dim, N> & morton, int level, map<Morton<Dim, N>, Cell>::const_iterator & begin, map<Morton<Dim, N>, Cell>::const_iterator & end) const
    {
        //std::cout << morton << " " << morton << std::endl;
        Morton<Dim, N> mortonBegin;
        Morton<Dim, N> mortonEnd;
        if (level >= N){
            begin= tree_.cbegin();
            end= tree_.cend();
            return;
        } else {
            mortonBegin = morton & masks_greater_equal_[level];
            mortonEnd = mortonBegin;
            mortonEnd.BitArrayBase<Morton<Dim,N>, Dim*N>::operator+=(masks_units_[level][0]);
        }
        //std::cout << mortonBegin << " " << mortonEnd << " " << level << std::endl;
        //std::cout << mortonBegin << " " << mortonEnd << " " << level << std::endl;
        begin = tree_.upper_bound(mortonBegin); --begin;
        end = tree_.lower_bound(mortonEnd);
        if (mortonEnd==0)
            end = tree_.end();
        if ( (begin->first & masks_greater_equal_[level]) != mortonBegin){
            end = tree_.end();
            begin = end;
        }
        //std::cout << begin->first << " " << end->first << std::endl;
    }

    inline Cell findCell(Morton<Dim, N> morton) const
    {
        auto itr = tree_.upper_bound(morton);
        --itr;
        Cell cell(itr->second);
        if ((morton & masks_greater_equal_[cell.level]) != itr->first){
            cell.level = -1;
            cell.ptr = NULL;
        }
        return cell;
    }

    map<Morton<Dim, N>, array<int, 2>> getFacesFine(int dir)
    {
        array<int, 2> level;
        map<Morton<Dim, N>, array<int, 2>> faces;
        auto itrF = faces.begin();
        Morton<Dim, N> morton, mortonNbr;
        Cell cellNbr;
        for (auto itr = tree_.cbegin(); itr != tree_.cend(); ++itr)
        {
            morton = itr->first;
            level[1] = itr->second.level;
            mortonNbr = morton - masks_units_[0][dir];
            cellNbr = findCell(mortonNbr);
            if ((cellNbr.level >= level[1]) || (cellNbr.level==-1))
            {
                level[0] = cellNbr.level;
                itrF = faces.insert(itrF, {morton, level});
            }

            level[0] = itr->second.level;
            mortonNbr = morton + masks_units_[level[0]][dir];
            cellNbr = findCell(mortonNbr);
            if ((cellNbr.level > level[0])|| (cellNbr.level==-1))
            {
                level[1] = cellNbr.level;
                itrF = faces.insert(itrF, {mortonNbr, level});
            }
        }
        return faces;
    }

    void balanceTree()
    {
/*
        std::map<Morton<Dim, N>, Cell> needRefinement;
        auto itrRef = needRefinement.begin();
        array<Morton<Dim, N>, 2> mortonNbr;
        array<Morton<Dim, N>, constexpr (1 << (Dim-1)) > mortonNbrSub;
        for (auto itr = tree_.cbegin(); itr != tree_.cend(); ++itr)
        {
            auto morton = itr->first;
            auto level = itr->second.level;
            for(auto i0=0; i0<Dim; ++i0)
            {
            mortonNbr[0] = morton - masks_units_[0][i0];
            mortonNbr[1] = morton + masks_units_[level][i0];
            for (auto j = 0; j < 2; ++j)
            {
                auto itrNbr = tree_.upper_bound(mortonNbr[j]);
                if (itrNbr == tree_.end())
                    continue;
                auto levelNbr = itrNbr->second.level;
                if ((mortonNbr[j] & masks_greater_than_[levelNbr]) != itrNbr->first) | (levelNbr >= level))
                    continue;
                if constexpr (Dim==1)
                {
                    if (level > (levelNbr+1))
                    {}                     
                }                
                else if constexpr (Dim==2)
                {

                }
                else if constexpr (Dim==3)
                {

                }
            }
        }
        return needRefinement;
        */
    }

    bool insert(const value_type &value);
    size_t deduceLevel(const Morton<Dim, N> &m) const
    {
        // logic: for a complete tree a cell is always part of an octant
        size_t level = m.countr_zero();
        if (level == 0)
            return level;
        size_t subLevel = level - 1;
        Morton<Dim, N> mSub = 1;
        mSub <<= subLevel;
        while (subLevel != -1)
        {
            if (tree_.count(m + mSub) == 0)
            {
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
            std::cout << (itr->first) << " " << itr->second << std::endl;
    }

    auto getCell(key_type morton, int level) const
    {
        std::array< std::array<int,Dim> ,2> cell;
        int dcell = (1<<level);
        for(auto i(0); i<Dim; ++i){
            cell[0][i] = encoder_.decode(morton,i).template getIntValue<int>();
            cell[1][i] = cell[0][i] + dcell;
        }
        return cell;
    }

    auto getCells() const{
        std::vector<std::array< std::array<int,Dim> ,2>> cells;
        cells.reserve(tree_.size());
        for(auto itr = tree_.cbegin(); itr != tree_.cend();++itr){
//            auto level = deduceLevel(itr->first);
            cells.push_back(getCell(itr->first, itr->second.level));
        }
        return cells;
    }
private:
    int levelCoarse_;
//    array<Morton<Dim, N>, N> masks_right_true_;
    array<Morton<Dim, N>, N> masks_greater_than_;
    array<Morton<Dim, N>, N> masks_greater_equal_;
    array<Morton<Dim, N>, N> masks_select_level_;
    array<Morton<Dim, N>, N> masks_unselect_level_;
    array<array<Morton<Dim, N>, Dim>, N> masks_units_;
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
//        masks_right_true_[i] = morton;
        masks_greater_than_[i] = ~morton;
        morton <<= 1;
    }
    morton = 0;
    morton = ~morton;
    for (auto i(0); i < N; ++i){
        masks_greater_equal_[i] = morton;
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


    morton = 1;
    for (auto i(0); i < N; ++i)
        for (auto j(0); j < Dim; ++j)
        {
            masks_units_[i][j] = morton;
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
bool Octree<Dim, N, T>::insert(const value_type &value)
{
    Morton<Dim, N> mortonNew(value.first);
    auto level = value.second.level;
    //std::cout << "mortonNew:    " << mortonNew << " " << level << std::endl;
    mortonNew &=  masks_greater_than_[level];
    Morton<Dim, N> mortonLoc(mortonNew);
    //std::cout << "mortonLoc:    " << mortonLoc << std::endl;
    mortonLoc.set(0, level); // put a one to be sure of the level
    //std::cout << "mortonLoc:    " << mortonLoc << std::endl;
    int levelParent = N;
    auto itr = tree_.begin();
    if (!tree_.empty())
    {
        itr = tree_.upper_bound(mortonLoc);
        if (itr != tree_.end())
            --itr;
        else {
            throw std::out_of_range("Cannot insert cell: parent not found.");
        }
        Morton<Dim, N> mortonParent= itr->first;
        levelParent= itr->second.level;
        if (levelParent <= level) return false;
        //std::cout << "mortonPar:    " << mortonParent << " " << levelParent << std::endl;

        if ((mortonNew & masks_greater_than_[levelParent-1]) != mortonParent){
                //std::cout << "masked:     " << (mortonNew & masks_greater_than_[levelParent]) << std::endl;
            throw std::out_of_range("Cannot insert cell: parent not found.");
        }

    }
    for (auto levelUp = levelParent-1; levelUp != (level-1) ; --levelUp){
        mortonLoc = mortonNew;
        for (auto i(0); i< numOctants; ++i){
            std::pair<Morton<Dim,N>, Cell> valueLoc;
            mortonLoc &= masks_greater_than_[levelUp];
            mortonLoc |= octants_[levelUp][i];
           // std::cout << "inserting: " << mortonLoc << std::endl;
           valueLoc.first = mortonLoc;
           valueLoc.second.level = levelUp;
           valueLoc.second.ptr = NULL;
            itr = tree_.insert(itr, valueLoc);
            //std::cout << "inserted2: " << mortonLoc << " " << levelUp << std::endl;
/*
            if (levelUp > levelCoarse_){
                auto mortonTest = mortonNew;
                mortonTest &= masks_greater_than_[4];
                auto itrTest = tree_.find(mortonTest);
                //std::cout << (itrTest == tree_.end()) << std::endl;
                //std::cout << mortonTest << std::endl;
                //std::cout << itrTest->first << std::endl;
                //std::cout << "Coarse cell :" << mortonLoc << " " << levelUp << std::endl;
            }
*/
        }
    }
    return true;
};

#endif // OCTREE_HPP