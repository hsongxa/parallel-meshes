/*
  Copyright (C) 2022 Hao Song

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef KD_TREE_H
#define KD_TREE_H

#include <cstddef> // std::size_t
#include <utility>
#include <memory>
#include <iterator>

namespace pmh {

struct kd_tree_node_base
{
   typedef kd_tree_node_base* base_ptr;

   kd_tree_node_base() : _parent(nullptr), _left(nullptr), _right(nullptr) {}

   base_ptr _parent;
   base_ptr _left;
   base_ptr _right;
};

template<typename P>
struct kd_tree_node : public kd_tree_node_base
{
   kd_tree_node(P p) : _value(p) {}
   P _value;
};

struct kd_tree_iterator_base
{
   typedef kd_tree_node_base::base_ptr     base_ptr;
   typedef std::size_t                     size_type;
   typedef std::bidirectional_iterator_tag iterator_category;
   typedef ptrdiff_t                       difference_type;

   kd_tree_iterator_base(base_ptr x, size_type d) : _node(x), _depth(d) {}

   // one important aspect of kd_tree is the splitting dimension associated
   // to each node - it is a function of the node's depth and is needed for
   // all the algorithms of kd_tree. we could store the depth (or splitting
   // dimension) in the node itself, which is a bit wasteful in memory, so
   // instead we store it here since most functions that return an iterator
   // can also generate the depth information and functions that use iterator
   // as an input also need to access the depth
   base_ptr  _node;
   size_type _depth;

   void increment()
   {
      // do not consider the special case of ++end(), which is
      // an undefined behavior according to the c++ standard 
      if (_node->_right != nullptr)
      {
         _node = _node->_right;
         ++_depth;
         while (_node->_left != nullptr)
         {
            _node = _node->_left;
            ++_depth;
         }
      }
      else
      {
         auto p = _node->_parent;
         while (_node == p->_right)
         {
            _node = p;
            --_depth;
            p = p->_parent;
         }
         _node = p;
         if (_depth > 0) --_depth;
      }
   }

   void decrement()
   {
      // do not consider the special case of --begin(), which is
      // an undefined behavior according to the c++ standard,
      // but --end() should give the last valid node when the tree
      // is not empty
      if (_node->_right == _node && _node->_left == _node && _node->_parent != nullptr)
      {
         // go to the rightmost
         _node = _node->_parent;
         while (_node->_right != nullptr)
         {
            _node = _node->_right;
            ++_depth;
         }
      }
      else if (_node->_left != nullptr)
      {
         _node = _node->_left;
         ++_depth;
         while (_node->_right != nullptr)
         {
            _node = _node->_right;
            ++_depth;
         }
      }
      else
      {
         auto p = _node->_parent;
         while (_node == p->_left)
         {
            _node = p;
            --_depth;
            p = p->_parent;
         }
         _node = p;
         if (_depth > 0) --_depth;
      }
   }
};

template<typename P, typename Ref, typename Ptr>
struct kd_tree_iterator : kd_tree_iterator_base
{
   typedef P   value_type;
   typedef Ref reference;
   typedef Ptr pointer;
   typedef kd_tree_iterator<P, Ref, Ptr> self;
   typedef kd_tree_node<P>* link_type;

   kd_tree_iterator(link_type x, size_type d) : kd_tree_iterator_base(x, d) {}
   kd_tree_iterator(const iterator& it) : kd_tree_iterator_base(it._node, it._depth) {}

   reference operator*() const { return static_cast<link_type>(_node)->_value; }
   pointer operator->() const { return &(operator*()); }

   self& operator++()
   {
      increment();
      return *this;
   }

   self operator++(int)
   {
      self tmp = *this;
      increment();
      return tmp;
   }

   self& operator--()
   {
      decrement();
      return *this;
   }

   self operator--(int)
   {
      self tmp = *this;
      decrement();
      return tmp;
   }
};

template<typename P, typename Ref, typename Ptr>
inline bool operator==(const kd_tree_iterator<P, Ref, Ptr>& it0, const kd_tree_iterator<P, Ref, Ptr>& it1)
{
   return it0._node == it1._node && it0._depth == it1._depth;
}

template<typename P, typename Ref, typename Ptr>
inline bool operator!=(const kd_tree_iterator<P, Ref, Ptr>& it0, const kd_tree_iterator<P, Ref, Ptr>& it1)
{
   return !(it0 == it1);
}

// TODO: implement KNN search
// TODO: implement orthogonal range search
// TODO: implement re-balancing
// reference algorithms: https://www.cs.umd.edu/class/fall2019/cmsc420-0201/Lects/lect13-kd-dict.pdf
//                       https://www.cs.umd.edu/class/fall2019/cmsc420-0201/Lects/lect14-kd-query.pdf
//                       https://www.geeksforgeeks.org/k-dimensional-tree/
// Note that
//    1. this implementation does NOT allow duplicate points;
//    2. the check for duplication is based on the implementation of Compare: two points are considered
//       duplicate if the Compare returns equal in ALL the dimensions from 0 to D - 1, and the user does
//       not need to provide a separate comparer for equality;
//    3. if a point's value in the splitting dimension is exactly the same as an existing node's splitting
//       value, the point is inserted to the RIGHT branch of that node, i.e., left children of a node are
//       STRICTLY smaller than it in the splitting dimension and the right children could be larger or
//       equal to the node's value in the splitting dimension.
// There exists a paper talking about handling duplicate points in k_d trees:
//    Meenakshi, Gill S. (2019) k-dLst Tree: k-d Tree with Linked List to Handle Duplicate Keys.
//    In: Rathore V., Worring M., Mishra D., Joshi A., Maheshwari S. (eds)
//        Emerging Trends in Expert Applicationsand Security.Advances in Intelligent Systemsand Computing,
//        vol 841. Springer, Singapore
template<std::size_t D, typename P, typename Compare, typename Alloc = std::allocator<P>>
class Kd_tree
{
public:
   typedef std::size_t                             size_type;
   typedef kd_tree_node<P>* link_type;
   typedef kd_tree_iterator<P, P&, P*>             iterator;
   typedef kd_tree_iterator<P, const P&, const P*> const_iterator;

private:
   typedef typename std::allocator_traits<Alloc>::template rebind_alloc<kd_tree_node<P>>  node_allocator;
   typedef typename std::allocator_traits<Alloc>::template rebind_traits<kd_tree_node<P>> node_alloc_traits;

   link_type get_node() { return node_alloc_traits::allocate(_allocator, 1); }

   void put_node(link_type node) { node_alloc_traits::deallocate(_allocator, node, 1); }

   link_type create_node(const P& p)
   {
      link_type tmp = get_node();
      try { node_alloc_traits::construct(_allocator, tmp, p); }
      catch (...) {
         put_node(tmp);
         throw;
      }
      return tmp;
   }

   void destroy_node(link_type node)
   {
      node_alloc_traits::destroy(_allocator, node); // destructor should not throw
      put_node(node);
   }

   link_type& root() const { return reinterpret_cast<link_type&>(_header->_parent); }

   iterator leftmost(link_type node, size_type depth) const
   {
      while (node->_left != nullptr)
      {
         node = static_cast<link_type>(node->_left);
         ++depth;
      }
      return iterator(node, depth);
   }

   void init()
   {
      _header = get_node();
      _header->_left = _header; // this configuration never changes
      _header->_right = _header; // this configuration never changes
      root() = nullptr;
   }

   static size_type dim_to_compare(size_type depth) { return depth % D; }

   static void custom_increment(iterator& it, size_type dim);

   bool are_equal(const P& p0, const P& p1) const
   {
      for (size_type i = 0; i < D; ++i)
         if (_compare(p0, p1, i) != 0)
            return false;
      return true;
   }

   iterator find_minimum(iterator start, size_type dim) const;

   std::pair<iterator, bool> insert_node(link_type node);

   void delete_node(iterator node);

private:
   node_allocator _allocator;
   Compare        _compare;
   size_type      _size;
   link_type      _header;

public:
   Kd_tree() : _allocator(), _compare(), _size(0), _header(nullptr) { init(); }
   ~Kd_tree() { clear(); put_node(_header); }

   static size_type dimension() { return D; }

   size_type size() const { return _size; }

   bool empty() const { return _size == 0; }

   iterator begin() const { return root() == nullptr ? end() : leftmost(root(), 0); }

   iterator end() const { return iterator(_header, 0); } // _header's depth is always zero (same as root)

   std::pair<iterator, bool> insert(const P& p)
   {
      link_type node = create_node(p);
      auto ret = insert_node(node);
      if (!ret.second) destroy_node(node);
      return ret;
   }

   void erase(iterator it) { return delete_node(it); }

   void erase(const P& p)
   {
      iterator it = find(p);
      if (it != end()) delete_node(it);
   }

   iterator find(const P& p) const;

   void clear();
};

template<std::size_t D, typename P, typename Compare, typename Alloc>
void Kd_tree<D, P, Compare, Alloc>::custom_increment(iterator& it, size_type dim)
{
   auto node = it._node;
   size_type depth = it._depth;
   if (node->_right != nullptr && dim_to_compare(depth) != dim)
   {
      node = node->_right;
      ++depth;
      while (node->_left != nullptr)
      {
         node = node->_left;
         ++depth;
      }
   }
   else
   {
      auto p = node->_parent;
      while (node == p->_right ||
         (node == p->_left && dim_to_compare(depth - 1) == dim))
      {
         node = p;
         --depth;
         p = p->_parent;
      }
      node = p;
      if (depth > 0) --depth;
   }
   it._node = node;
   it._depth = depth;
}

template<std::size_t D, typename P, typename Compare, typename Alloc>
typename Kd_tree<D, P, Compare, Alloc>::iterator
Kd_tree<D, P, Compare, Alloc>::find_minimum(iterator start, size_type dim) const
{
   iterator it = leftmost(static_cast<link_type>(start._node), start._depth);
   iterator min = it;
   for (custom_increment(it, dim); it._depth >= start._depth; custom_increment(it, dim))
      if (_compare(*it, *min, dim) < 0) min = it;
   return min;
}

template<std::size_t D, typename P, typename Compare, typename Alloc>
std::pair<typename Kd_tree<D, P, Compare, Alloc>::iterator, bool>
Kd_tree<D, P, Compare, Alloc>::insert_node(link_type node)
{
   link_type curr = root();
   if (curr == nullptr)
   {
      root() = node;
      node->_parent = _header;
      ++_size;
      return std::make_pair(iterator(node, 0), true);
   }
   else
   {
      size_type depth = 0;
      while (true)
      {
         if (are_equal(node->_value, curr->_value))
            return std::make_pair(iterator(curr, depth), false); // duplicate

         if (_compare(node->_value, curr->_value, dim_to_compare(depth)) < 0)
         {
            if (curr->_left == nullptr)
            {
               curr->_left = node;
               node->_parent = curr;
               ++depth;
               ++_size;
               break;
            }
            else
               curr = static_cast<link_type>(curr->_left);
         }
         else
         {
            if (curr->_right == nullptr)
            {
               curr->_right = node;
               node->_parent = curr;
               ++depth;
               ++_size;
               break;
            }
            else
               curr = static_cast<link_type>(curr->_right);
         }
         ++depth;
      }
      return std::make_pair(iterator(node, depth), true);
   }
}

template<std::size_t D, typename P, typename Compare, typename Alloc>
void Kd_tree<D, P, Compare, Alloc>::delete_node(iterator node)
{
   auto curr = static_cast<link_type>(node._node);
   size_type depth = node._depth;
   while (curr != nullptr)
   {
      if (curr->_right != nullptr) // replace curr by the min of its right branch
      {
         auto min = find_minimum(iterator(static_cast<link_type>(curr->_right), depth + 1), dim_to_compare(depth));
         curr->_value = *min; // assignment operation of P
         curr = static_cast<link_type>(min._node);
         depth = min._depth;
      }
      else if (curr->_left != nullptr) // replace curr by the min of its left branch and flip (left => right)
      {
         auto min = find_minimum(iterator(static_cast<link_type>(curr->_left), depth + 1), dim_to_compare(depth));
         curr->_value = *min; // assignment operation of P
         curr->_right = curr->_left;
         curr->_left = nullptr;
         curr = static_cast<link_type>(min._node);
         depth = min._depth;
      }
      else // delete - this is a leaf node
      {
         if (curr == root())
            root() = nullptr;
         else if (curr->_parent->_left == curr)
            curr->_parent->_left = nullptr;
         else
            curr->_parent->_right = nullptr;

         destroy_node(curr);
         --_size;

         curr = nullptr;
      }
   }
}

template<std::size_t D, typename P, typename Compare, typename Alloc>
typename Kd_tree<D, P, Compare, Alloc>::iterator
Kd_tree<D, P, Compare, Alloc>::find(const P& p) const
{
   link_type curr = root();
   size_type depth = 0;
   while (curr != nullptr)
   {
      if (are_equal(curr->_value, p))
         return iterator(curr, depth);

      if (_compare(p, curr->_value, dim_to_compare(depth)) < 0)
         curr = static_cast<link_type>(curr->_left);
      else
         curr = static_cast<link_type>(curr->_right);
      ++depth;
   }
   return end();
}

template<std::size_t D, typename P, typename Compare, typename Alloc>
void Kd_tree<D, P, Compare, Alloc>::clear()
{
   link_type curr = root();
   while (curr != nullptr && curr != _header)
   {
      if (curr->_left != nullptr)
         curr = static_cast<link_type>(curr->_left);
      else if (curr->_right != nullptr)
         curr = static_cast<link_type>(curr->_right);
      else
      {
         if (curr == root())
            root() = nullptr;
         else if (curr->_parent->_left == curr)
            curr->_parent->_left = nullptr;
         else
            curr->_parent->_right = nullptr;

         auto parent = static_cast<link_type>(curr->_parent);
         destroy_node(curr);
         --_size;
         curr = parent;
      }
   }
}

} // namespace pmh

#endif
