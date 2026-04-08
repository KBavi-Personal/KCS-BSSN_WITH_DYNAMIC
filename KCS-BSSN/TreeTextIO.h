#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>

#include "Tree.h"
#include "TreeNode.h"
#include "SampleData.h"
#include "BipartiteNetwork.h"

class TreeTextIO
{
public:
    
    static bool Save(const Tree* tree, const std::string& path)
    {
        if (!tree || !tree->sd || !tree->root) return false;

        std::ofstream out(path);
        if (!out.is_open()) return false;

        out << "TREE_V3\n";


        out << "NODE_COUNT " << tree->all_nodes.size() << "\n";
        out << "ROOT_ID " << tree->root->id << "\n";
        out << "NUM_OF_NODES " << tree->num_of_nodes << "\n";


        out << "P_INDEX_VEC " << tree->P_index_vec_for_nodes.size();
        for (User* u : tree->P_index_vec_for_nodes)
            out << " " << (u ? u->id : -1);
        out << "\n";


        out << "NODES_LIST " << tree->Nodes.size();
        for (TreeNode* n : tree->Nodes)
            out << " " << (n ? n->id : -1);
        out << "\n";

        out << "ALL_NODES_LIST " << tree->all_nodes.size();
        for (TreeNode* n : tree->all_nodes)
            out << " " << (n ? n->id : -1);
        out << "\n";

        out << "USER_LEAF " << tree->user_leaf.size();
        for (TreeNode* leaf : tree->user_leaf)
            out << " " << (leaf ? leaf->id : -1);
        out << "\n";


        for (TreeNode* n : tree->all_nodes)
        {
            if (!n) return false;

            out << "NODE_BEGIN\n";

            out << "ID " << n->id << "\n";
            out << "DEPTH " << n->depth << "\n";
            out << "PIVOT_USER_INDEX " << n->myUserIndexForPivot << "\n";
            out << "IS_LEAF " << (n->is_leaf ? 1 : 0) << "\n";
            out << "IS_ROOT " << (n->isroot ? 1 : 0) << "\n";
            out << "TOTAL_N_KEYS " << n->total_n_keys << "\n";
            out << "UBS " << n->ub_sub << " " << n->ub_w_in << " " << n->ub_w_out << "\n";
            out << "P_S_P_R " << n->P_s_size << " " << n->P_r_size << "\n";
            out << "PARENT_ID " << (n->parent ? n->parent->id : -1) << "\n";


            out << "DIST_PS " << (2 * n->P_s_size);
            for (int i = 0; i < 2 * n->P_s_size; ++i)
                out << " " << n->min_max_dist_P_s[i];
            out << "\n";


            out << "USERS " << n->users.size();
            for (User* u : n->users)
                out << " " << (u ? u->id : -1);
            out << "\n";

            out << "CHILDREN " << n->children.size();
            for (TreeNode* c : n->children)
                out << " " << (c ? c->id : -1);
            out << "\n";


            out << "KEYS_UB_FSUM " << n->Keys_ub_Fsum.size();
            for (double x : n->Keys_ub_Fsum)
                out << " " << x;
            out << "\n";


            out << "KEYS_UB_FMAX " << n->Keys_ub_Fmax.size();
            for (double x : n->Keys_ub_Fmax)
                out << " " << x;
            out << "\n";

            out << "NODE_END\n";
        }

        out << "TREE_END\n";
        return true;
    }


    static Tree* Load(SampleData* sd, const std::string& path)
    {
        if (!sd || !sd->bn) return nullptr;

        std::ifstream in(path);
        if (!in.is_open()) {
            std::cout << "Could not open file: " << path << "\n";
            return nullptr;
        }

        std::string tok;
        in >> tok;
        if (tok != "TREE_V3") {
            std::cout << "Unsupported tree file version.\n";
            return nullptr;
        }

        size_t node_count = 0;
        int root_id = -1;
        int num_of_nodes = 0;


        if (!(in >> tok) || tok != "NODE_COUNT") return nullptr;
        in >> node_count;

        if (!(in >> tok) || tok != "ROOT_ID") return nullptr;
        in >> root_id;

        if (!(in >> tok) || tok != "NUM_OF_NODES") return nullptr;
        in >> num_of_nodes;

        Tree* tree = new Tree(sd);
        tree->inx = nullptr;
        tree->num_of_nodes = num_of_nodes;


        if (!(in >> tok) || tok != "P_INDEX_VEC") { delete tree; return nullptr; }
        size_t pivN = 0;
        in >> pivN;
        tree->P_index_vec_for_nodes.clear();
        tree->P_index_vec_for_nodes.reserve(pivN);
        for (size_t i = 0; i < pivN; ++i) {
            int uid;
            in >> uid;
            if (uid >= 0 && uid < (int)sd->bn->users_vec.size())
                tree->P_index_vec_for_nodes.push_back(sd->bn->users_vec[uid]);
            else
                tree->P_index_vec_for_nodes.push_back(nullptr);
        }


        if (!(in >> tok) || tok != "NODES_LIST") { delete tree; return nullptr; }
        size_t nodesListN = 0;
        in >> nodesListN;
        std::vector<int> nodes_list_ids(nodesListN);
        for (size_t i = 0; i < nodesListN; ++i) in >> nodes_list_ids[i];


        if (!(in >> tok) || tok != "ALL_NODES_LIST") { delete tree; return nullptr; }
        size_t allListN = 0;
        in >> allListN;
        std::vector<int> all_nodes_ids(allListN);
        for (size_t i = 0; i < allListN; ++i) in >> all_nodes_ids[i];


        if (!(in >> tok) || tok != "USER_LEAF") { delete tree; return nullptr; }
        size_t userLeafN = 0;
        in >> userLeafN;
        std::vector<int> user_leaf_ids(userLeafN);
        for (size_t i = 0; i < userLeafN; ++i) in >> user_leaf_ids[i];


        tree->all_nodes.clear();
        tree->all_nodes.reserve(node_count);

        std::unordered_map<int, TreeNode*> by_id;
        by_id.reserve(node_count * 2);

        std::vector<int> read_ids;
        read_ids.reserve(node_count);

        std::vector<std::vector<int>> children_ids(node_count);
        std::vector<std::vector<int>> users_ids(node_count);
        std::vector<int> parent_ids(node_count, -1);


        for (size_t idx = 0; idx < node_count; ++idx)
        {
            if (!(in >> tok) || tok != "NODE_BEGIN") {
                delete tree;
                return nullptr;
            }

            int id = -1;
            int depth = 0;
            int pivot_idx = -1;
            int is_leaf_int = 0;
            int is_root_int = 0;
            int total_n_keys = 0;
            float ub_sub = 0.0f, ub_w_in = 0.0f, ub_w_out = 0.0f;
            int P_s_size = 0, P_r_size = 0;
            int parent_id = -1;

            if (!(in >> tok) || tok != "ID") { delete tree; return nullptr; }
            in >> id;

            if (!(in >> tok) || tok != "DEPTH") { delete tree; return nullptr; }
            in >> depth;

            if (!(in >> tok) || tok != "PIVOT_USER_INDEX") { delete tree; return nullptr; }
            in >> pivot_idx;

            if (!(in >> tok) || tok != "IS_LEAF") { delete tree; return nullptr; }
            in >> is_leaf_int;

            if (!(in >> tok) || tok != "IS_ROOT") { delete tree; return nullptr; }
            in >> is_root_int;

            if (!(in >> tok) || tok != "TOTAL_N_KEYS") { delete tree; return nullptr; }
            in >> total_n_keys;

            if (!(in >> tok) || tok != "UBS") { delete tree; return nullptr; }
            in >> ub_sub >> ub_w_in >> ub_w_out;

            if (!(in >> tok) || tok != "P_S_P_R") { delete tree; return nullptr; }
            in >> P_s_size >> P_r_size;

            if (!(in >> tok) || tok != "PARENT_ID") { delete tree; return nullptr; }
            in >> parent_id;

            TreeNode* n = new TreeNode(depth, P_s_size, P_r_size, tree, (is_leaf_int != 0), pivot_idx, total_n_keys);
            n->id = id;
            n->depth = depth;
            n->isroot = (is_root_int != 0);
            n->ub_sub = ub_sub;
            n->ub_w_in = ub_w_in;
            n->ub_w_out = ub_w_out;


            if (!(in >> tok) || tok != "DIST_PS") { delete tree; return nullptr; }
            int distCount = 0;
            in >> distCount;
            if (distCount != 2 * P_s_size) {
                delete n;
                delete tree;
                return nullptr;
            }
            for (int i = 0; i < distCount; ++i)
                in >> n->min_max_dist_P_s[i];


            if (!(in >> tok) || tok != "USERS") { delete tree; return nullptr; }
            size_t uN = 0;
            in >> uN;
            users_ids[idx].resize(uN);
            for (size_t i = 0; i < uN; ++i)
                in >> users_ids[idx][i];

            if (!(in >> tok) || tok != "CHILDREN") { delete tree; return nullptr; }
            size_t cN = 0;
            in >> cN;
            children_ids[idx].resize(cN);
            for (size_t i = 0; i < cN; ++i)
                in >> children_ids[idx][i];


            if (!(in >> tok) || tok != "KEYS_UB_FSUM") { delete tree; return nullptr; }
            size_t fsumN = 0;
            in >> fsumN;
            n->Keys_ub_Fsum.resize(fsumN);
            for (size_t i = 0; i < fsumN; ++i)
                in >> n->Keys_ub_Fsum[i];


            if (!(in >> tok) || tok != "KEYS_UB_FMAX") { delete tree; return nullptr; }
            size_t fmaxN = 0;
            in >> fmaxN;
            n->Keys_ub_Fmax.resize(fmaxN);
            for (size_t i = 0; i < fmaxN; ++i)
                in >> n->Keys_ub_Fmax[i];

            if (!(in >> tok) || tok != "NODE_END") {
                delete tree;
                return nullptr;
            }

            tree->all_nodes.push_back(n);
            by_id[id] = n;
            read_ids.push_back(id);
            parent_ids[idx] = parent_id;
        }

        if (!(in >> tok) || tok != "TREE_END") {
            delete tree;
            return nullptr;
        }


        for (size_t idx = 0; idx < node_count; ++idx)
        {
            TreeNode* n = by_id[read_ids[idx]];

            n->children.clear();
            n->children.reserve(children_ids[idx].size());
            for (int cid : children_ids[idx]) {
                auto it = by_id.find(cid);
                if (it != by_id.end())
                    n->children.push_back(it->second);
            }

            n->users.clear();
            n->users.reserve(users_ids[idx].size());
            for (int uid : users_ids[idx]) {
                if (uid >= 0 && uid < (int)sd->bn->users_vec.size())
                    n->users.push_back(sd->bn->users_vec[uid]);
                else
                    n->users.push_back(nullptr);
            }

            if (parent_ids[idx] >= 0) {
                auto pit = by_id.find(parent_ids[idx]);
                n->parent = (pit != by_id.end()) ? pit->second : nullptr;
            }
            else {
                n->parent = nullptr;
            }
        }


        tree->Nodes.clear();
        tree->Nodes.reserve(nodes_list_ids.size());
        for (int nid : nodes_list_ids) {
            auto it = by_id.find(nid);
            if (it != by_id.end())
                tree->Nodes.push_back(it->second);
        }


        auto rit = by_id.find(root_id);
        tree->root = (rit != by_id.end()) ? rit->second : nullptr;
        if (!tree->root) {
            delete tree;
            return nullptr;
        }


        if (all_nodes_ids.size() == tree->all_nodes.size()) {
            std::vector<TreeNode*> reordered;
            reordered.reserve(all_nodes_ids.size());
            for (int nid : all_nodes_ids) {
                auto it = by_id.find(nid);
                if (it != by_id.end())
                    reordered.push_back(it->second);
            }
            if (reordered.size() == tree->all_nodes.size())
                tree->all_nodes.swap(reordered);
        }

        tree->user_leaf.clear();
        tree->user_leaf.resize(userLeafN, nullptr);
        for (size_t i = 0; i < userLeafN; ++i) {
            int nid = user_leaf_ids[i];
            if (nid >= 0) {
                auto it = by_id.find(nid);
                if (it != by_id.end())
                    tree->user_leaf[i] = it->second;
            }
        }

        return tree;
    }
};