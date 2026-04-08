#pragma once
#include <vector>
#include <map>
#include <utility>
#include "User.h"

class BipartiteNetwork
{
public:
    int no_of_users;
    std::vector<User*> users_vec; // users of social network

    BipartiteNetwork();
    ~BipartiteNetwork();

    void Reset_Users_pruned_flag()
    {
        for (auto u : users_vec)
            u->isPruned = false;
    }

    // ---------------- Edge support API ----------------
    // Always stored as (u>=v)
    int  Get_edge_sup(int u, int v) const;
    void Set_edge_sup(int u, int v, int sup);
    void Inc_edge_sup(int u, int v, int delta = 1);
    bool Has_edge(int u, int v) const;

    void Add_edge_if_missing(int u, int v);

    void Clear_edges_sup();
    void Reset_edges_sup_values();
    const std::map<std::pair<int, int>, int>& Edges_sup() const;
    std::map<std::pair<int, int>, int>& Edges_sup_mut();

    std::string print();

private:
    std::map<std::pair<int, int>, int> edges_sup; // e(u,v)=sup(e)

    static inline std::pair<int, int> Norm_edge(int u, int v)
    {
        return (u >= v) ? std::make_pair(u, v) : std::make_pair(v, u);
    }
};
