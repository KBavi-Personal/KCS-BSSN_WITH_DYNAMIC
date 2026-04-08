#include "BipartiteNetwork.h"
#include "User.h"

BipartiteNetwork::BipartiteNetwork()
{
}
BipartiteNetwork::~BipartiteNetwork() {
	for (int i = 0; i < users_vec.size(); i++)
		delete users_vec[i];
}
std::string BipartiteNetwork::print() {
	std::string text = "Bipartite Network:\n";
	//printf("========== BipartiteNetwork =======================================\n");
	int graph_size = users_vec.size();
	text += "Social Network... # of users=" + std::to_string(graph_size) + "  \n";
	for (int i = 0; i < graph_size; i++) {
		text += "User[" + std::to_string(i) + "]:  id=" + std::to_string(users_vec[i]->id) + "  #inNeighbors=" + std::to_string(users_vec[i]->inNeighbors.size()) + "  #outNeighbors=" + std::to_string(users_vec[i]->outNeighbors.size()) + "   #checkin_loc=" + std::to_string(users_vec[i]->checkin_locations.size()) + "\n";
		text += "Checkin Loc - frequency ";
		for (auto& u : users_vec[i]->checkin_locations) {
			text += "(" + std::to_string(u.first) + " - " + std::to_string(u.second) + ") ";
		}
		text += "\ninNeighbors - influnce";
		for (auto& u : users_vec[i]->inNeighbors) {
			text += "(" + std::to_string(u.first) + " - " + std::to_string(u.second) + ")  ";
		}
		text += "\noutNeighbors - influnce";
		for (auto& u : users_vec[i]->outNeighbors) {
			text += "(" + std::to_string(u.first) + " - " + std::to_string(u.second) + ")  ";
		}
		text += "\n\n";
	}


	text += "outNeighbors\n";
	for (int i = 0; i < graph_size; i++) {
		for (auto out : users_vec[i]->outNeighbors)
			text += std::to_string(users_vec[i]->id) + "  " + std::to_string(out.first) + " " + std::to_string(out.second) + "\n";
	}
	text += "inNeighbors\n";
	for (int i = 0; i < graph_size; i++) {
		for (auto in : users_vec[i]->inNeighbors)
			text += std::to_string(users_vec[i]->id) + "  " + std::to_string(in.first) + " " + std::to_string(in.second) + "\n";

	}
	text += "Edges and sup\n";
	for (const auto& u : edges_sup)
		text += "sup(" + std::to_string(u.first.first) + "," + std::to_string(u.first.second) + ")=" + std::to_string(u.second) + " \n";

	return text;
}

int BipartiteNetwork::Get_edge_sup(int u, int v) const
{
	auto it = edges_sup.find(Norm_edge(u, v));
	return (it == edges_sup.end()) ? 0 : it->second;
}

void BipartiteNetwork::Set_edge_sup(int u, int v, int sup)
{
	edges_sup[Norm_edge(u, v)] = sup;
}

void BipartiteNetwork::Inc_edge_sup(int u, int v, int delta)
{
	edges_sup[Norm_edge(u, v)] += delta;
}

bool BipartiteNetwork::Has_edge(int u, int v) const
{
	return edges_sup.count(Norm_edge(u, v)) == 1;
}

void BipartiteNetwork::Add_edge_if_missing(int u, int v)
{
	auto key = Norm_edge(u, v);
	if (edges_sup.find(key) == edges_sup.end())
		edges_sup.emplace(key, 0);
}

void BipartiteNetwork::Clear_edges_sup()
{
	edges_sup.clear();
}

void BipartiteNetwork::Reset_edges_sup_values()
{
	for (auto& e : edges_sup) e.second = 0;
}

const std::map<std::pair<int, int>, int>& BipartiteNetwork::Edges_sup() const
{
	return edges_sup;
}

std::map<std::pair<int, int>, int>& BipartiteNetwork::Edges_sup_mut()
{
	return edges_sup;
}


