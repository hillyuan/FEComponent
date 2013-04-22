/*==============================================================================

                                    FEComponent

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 2013 Xi YUAN

   This file is part of FECompnent.

   FEComponent is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FEComponent is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with FEComponent. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                        Data structure

  ==============================================================================*/

#ifndef __VECTOR_TREE_H
#define __VECTOR_TREE_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace XYLIB
{

class CVectorTree
{
    typedef vector<double_t>               ValueType;
    typedef vector<double_t>               Row;

public:
    struct Node {
        typedef vector<Node*>                       Children;
        typedef typename Children::iterator         Iterator;
        typedef typename Children::const_iterator   const_Iterator;

        ValueType  _value;
        Children   _children;

        Node(const ValueType value): _value(value) {}
    };

    typedef typename Node::Children         Children;
    typedef typename Node::Iterator         Iterator;
    typedef typename Node::const_Iterator   ConstIterator;

public:
    explicit CVectorTree(const ValueType& value) {
        m_ncomp = 1;
        m_root = new Node(value);
    }

    ~CVectorTree() {
        Clear();
        delete m_root;
    }

    Node*       GetRoot() {return m_root;}
    const Node* GetRoot() const {return m_root;}

    static Node* Append(Node* node, const ValueType& value) {
        Node* child = new Node(value);
        node->_children.push_back(child);
        return child;
    }


    // Delete all nodes except root node
    void Clear() {
        Clear(m_root);
    }

    void Show(const Node* node, int depth=0) {
        unsigned int i;
        cout << "Level:" << depth << std::endl;
        if( !node->_value.empty() ) {
            for(i=0;i<node->_value.size();i++) {
                std::cout << i << ":" << node->_value[i] << std::endl;
            }
        }
        const Children& children = node->_children;
        for(ConstIterator it=children.begin(); it!=children.end(); ++it) {
            Show(*it,depth+1);
        }
    }

    bool isConsistent() const {
        return isConsistent(m_root);
    }

    bool GetValue(ValueType& idata,ValueType& odata ) {
        int clevel = idata.size();
        if( idata.empty() ) {
            return GetValue(m_root, clevel, odata);
        }
        return GetValue(m_root, idata, clevel, odata);
    }

    bool GetGrad(ValueType& idata,ValueType& odata) {
        int clevel = idata.size();
        if( idata.empty() ) {
            return GetValue(m_root, clevel, odata);
        }
        return GetGrad(m_root, idata, clevel, odata);
    }

    CVectorTree(int ncomp, std::ifstream& file) {
        std::vector<Row> table;

        while(!file.eof()){
            std::string line;
            std::getline(file, line);
            std::istringstream is(line);
            Row row;
            while(!is.eof()) {
                double_t data;
                is >> data;
                row.push_back(data);
            }
            table.push_back(row);
        }
        m_ncomp = ncomp;

        Node *node;
        unsigned int ncol=table[0].size();
        unsigned int maxcol=ncol;
        std::vector<unsigned int> index;
        index.push_back(0);
        index.push_back(table.size());
        RecursiveParsing(node, table, ncol, index[0], index[1], maxcol);
    }

private:
    void Clear(Node* node) {
        Children& children = node->_children;
        for(Iterator it=children.begin(); it!=children.end(); ++it) {
            Clear(*it);
        }
        children.clear();

        if(node!=m_root) {delete node;}
    }

    bool isConsistent(Node* node) const {
        if( node->_value.size() != node->_children.size()) return false;
        Children& children = node->_children;
        for(Iterator it=children.begin(); it!=children.end(); ++it) {
            if( isConsistent(*it)==false ) return false;
        }
        return true;
    }

    bool GetValue(Node* node, ValueType& idata, int& clevel, ValueType& odata ) {
        int level1;
        unsigned int k;
        double_t lambda;
        vector<double_t>::iterator low, low1;
        ValueType odata1;

        low = find(node->_value.begin(), node->_value.end(), idata[clevel-1]);
        if( low!=node->_value.end() ) {
            cout << "Found:" << clevel << "  "<< *low << endl;
            level1 = low- node->_value.begin();
            if(clevel==1) {
                for(k=0; k<m_ncomp; k++) {
                    odata.push_back(node->_children[level1]->_value[k]);
                }
                return true;
            } else {
                clevel--;
                Node* child = node->_children[level1];
                if( !child ) return false;
                GetValue(child, idata, clevel, odata);
            }
        } else {
            low=std::lower_bound(node->_value.begin(), node->_value.end(), idata[clevel-1]);
            level1 = low- node->_value.begin();
            cout << "Found low:" << clevel << "  "<< idata[clevel-1] << "  " << level1 << endl;
            if(clevel==1) {
                if( low==node->_value.begin()) {
                    for(k=0; k<m_ncomp; k++) {
                        odata.push_back(node->_children[level1]->_value[k]);
                    }
                } else if ( low==node->_value.end()) {  // not found
                    for(k=0; k<m_ncomp; k++) {
                        odata.push_back(node->_children[level1-1]->_value[k]);
                    }
                } else {
                    low1 = low-1;
                    lambda = (*low-idata[clevel-1])/(*low-*low1);
                    cout << "lambda=" << lambda << "  " << *low1 << "  " << *low << "   " << idata[clevel-1] << endl;
                    for(k=0; k<m_ncomp; k++) {
                        cout << node->_children[level1-1]->_value[k] << "  " << node->_children[level1]->_value[k] << endl;
                        odata.push_back(node->_children[level1-1]->_value[k] + lambda*
                            (node->_children[level1]->_value[k]-node->_children[level1-1]->_value[k]) );
                    }
                }
                return true;
            } else {
                low1 = low-1;
                lambda = (*low-idata[clevel-1])/(*low-*low1);
                cout << "lambda=" << lambda << "  " << *low1 << "  " << *low << "   " << idata[clevel-1] << endl;
                clevel--;
                Node* child = node->_children[level1-1];
                if( !child ) return false;
                GetValue(child, idata, clevel, odata);
                Node* child1 = node->_children[level1];
                if( !child1 ) return false;
                GetValue(child1, idata, clevel, odata1);
                for(k=0; k<m_ncomp; k++) {
                    odata[k] = odata[k] + lambda* (odata1[k]-odata[k]);
                }
            }
        }
        return true;
    }

    bool GetValue(Node* node, int& clevel, ValueType& odata ) {
        clevel--;
        if(clevel==1) {
            for(unsigned int k=0; k<m_ncomp; k++) {
                odata.push_back(node->_children[0]->_value[k]);
            }
            return true;
        }
        Children& children = node->_children;
        GetValue(children[0], clevel, odata );

        return  true;
    }

    bool GetGrad(Node* node, ValueType& idata, int& clevel, ValueType& odata ) {
        int level1;
        unsigned int k;
        double_t a1,a2,b1,b2;
        vector<double_t>::iterator low, low1;
        ValueType odata1;

        low = find(node->_value.begin(), node->_value.end(), idata[clevel-1]);
        if( low!=node->_value.end() ) {
            cout << "Grad found:" << clevel << "  "<< *low << endl;
            level1 = low- node->_value.begin();
            if(clevel==1) {
                if( low==node->_value.begin() ) {
                    a1 = *low;
                    low1 = low+1; a2=*low1;
                    if( a2<=a1 ) {
                        cout << "Math error in GetGrad" << endl;
                        return false;
                    }
                    for(k=0; k<m_ncomp; k++) {
                        b1 = node->_children[level1]->_value[k];
                        b2 = node->_children[level1+1]->_value[k];
                        odata.push_back((b2-b1)/(a2-a1));
                    }
                }  else {
                    a2 = *low;
                    low1 = low-1; a1=*low1;
                    if( a2<=a1 ) {
                        cout << "Math error in GetGrad" << endl;
                        return -1;
                    }
                    for(k=0; k<m_ncomp; k++) {
                        b1 = node->_children[level1-1]->_value[k];
                        b2 = node->_children[level1]->_value[k];
                        odata.push_back((b2-b1)/(a2-a1));
                    }
                }
                return true;
            } else {
                clevel--;
                Node* child = node->_children[level1];
                if( !child ) return false;
                GetGrad(child, idata, clevel, odata);
            }
        } else {
            low=std::lower_bound(node->_value.begin(), node->_value.end(), idata[clevel-1]);
            level1 = low- node->_value.begin();
            cout << "Grad found low:" << clevel << "  "<< idata[clevel-1] << "  " << level1 << endl;
            if(clevel==1) {
                if( low==node->_value.begin()) {
                /*    a1 = *low;
                    low1 = low+1; a2=*low1;
                    if( a2<=a1 ) {
                        cout << "Math error in GetGrad" << endl;
                        return;
                    }*/
                    for(k=0; k<m_ncomp; k++) {
                   /*     b1 = node->_children[level1]->_value[k];
                        b2 = node->_children[level1+1]->_value[k];
                        odata.push_back((b2-b1)/(a2-a1));*/
                        odata.push_back(0.0);
                    }
                } else if ( low==node->_value.end()) {  // not found
                  /*  low1 = low-1; a1=*low1;
                    low1=low1-1;  a2 =*low1;
                    if( a2<=a1 ) {
                        cout << "Math error in GetGrad" << endl;
                        return;
                    }*/
                    for(k=0; k<m_ncomp; k++) {
                     /*   b1 = node->_children[level1-2]->_value[k];
                        b2 = node->_children[level1-1]->_value[k];
                        odata.push_back((b2-b1)/(a2-a1));*/
                        odata.push_back(0.0);
                    }
                } else {
                    a2 = *low;
                    low1 = low-1; a1=*low1;
                    if( a2<=a1 ) {
                        cout << "Math error in GetGrad" << endl;
                        return -1;
                    }
                    for(k=0; k<m_ncomp; k++) {
                        b1 = node->_children[level1-1]->_value[k];
                        b2 = node->_children[level1]->_value[k];
                        odata.push_back((b2-b1)/(a2-a1));
                    }
                }
                return true;
            } else {
                low1 = low-1;
                double_t lambda = (*low-idata[clevel-1])/(*low-*low1);
                cout << "lambda=" << lambda << "  " << *low1 << "  " << *low << "   " << idata[clevel-1] << endl;
                clevel--;
                Node* child = node->_children[level1-1];
                if( !child ) return false;
                GetGrad(child, idata, clevel, odata);
                Node* child1 = node->_children[level1];
                if( !child1 ) return false;
                GetGrad(child1, idata, clevel, odata1);
                for(k=0; k<m_ncomp; k++) {
                    odata[k] = odata[k] + lambda* (odata1[k]-odata[k]);
                }
            }
        }
        return true;
    }

    void RecursiveParsing(Node* node, std::vector<Row>& table, unsigned int& ncol, unsigned int& sindex, unsigned int& eindex,
                          unsigned int& maxcol) {
        unsigned int i,j,k,ncol1;
        ValueType data;
        vector<unsigned int> index1;

        cout << "Parsing: col="  << ncol << "  row=" << sindex<< ", "<< eindex << std::endl;

        index1.push_back(sindex);

        if( ncol==m_ncomp ) { // leaf
            for(k=0; k<m_ncomp; k++) {
                data.push_back(table[sindex][ncol-m_ncomp+k]);
            }
            node=Append(node, data);
            return;
        } else {
            data.push_back(table[sindex][ncol-1]);
            for(j=sindex+1; j<eindex; j++) {
                if(ncol==m_ncomp) {

                } else if(table[j][ncol-1]!=table[j-1][ncol-1]) {
                    data.push_back(table[j][ncol-1]);
                    index1.push_back(j);
                }
            }
        }
        if(ncol==maxcol) {
                m_root = new Node(data);
                node = m_root;
        } else {
                node=Append(node, data);
        }


        index1.push_back(eindex);
        if( ncol==m_ncomp ) return;

        /* Parsing next level */
        ncol1=ncol-1;
        for(i=0; i<index1.size()-1; i++) {
            RecursiveParsing(node, table, ncol1, index1[i], index1[i+1], maxcol);
        }
    }

private:
    Node* m_root;
    unsigned int m_ncomp;
};

}

#endif
