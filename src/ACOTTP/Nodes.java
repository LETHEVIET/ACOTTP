package de.adrianwilke.acotspjava;

import java.util.ArrayList;

class Node_Base<NodeT>
{
    public NodeT parent = null;
    public ArrayList<NodeT> childs = new ArrayList<NodeT>(2);
}

class Node extends Node_Base<Node> {

    public double _past_pheromone; // pheromone of the edge from the parent node to itself at the time of last local evaporation
    public double _heuristic;      // heuristic of the edge from the parent node to itself
    public int _local_restart_times;
    public int _local_evap_times;
    public int _n_child_leaf;

public Node(Node parent, double heuristic, int n_child_leaf)
{
    this.parent = parent;
    _local_evap_times = 0;
    _local_restart_times = 0;
    _heuristic = heuristic;
    _n_child_leaf = n_child_leaf;
}

public void local_restart_if_needed(int global_restart_times, double past_trail_restart)
{
    assert(_local_restart_times <= global_restart_times);
    if (_local_restart_times < global_restart_times)
    {
        _past_pheromone = past_trail_restart;
        _local_restart_times = global_restart_times;
        _local_evap_times = 0;
    }
}

    public double get_pheromone(double one_minus_rho,
                           double past_trail_restart,
                           double past_trail_min,
                           int global_restart_times, int global_evap_times)
{
    double pheromone;

    local_restart_if_needed(global_restart_times, past_trail_restart);
    assert(_local_evap_times <= global_evap_times);
    pheromone = _past_pheromone;
    assert(pheromone > 0);

    if (_local_evap_times < global_evap_times)
    {
        pheromone *= Math.pow(one_minus_rho, global_evap_times - _local_evap_times);
        if (pheromone < past_trail_min)
            pheromone = past_trail_min;
    }

    return pheromone;
}

    public void pay_evaporation_debt(double one_minus_rho,
                                double past_trail_restart,
                                double past_trail_min,
                                int global_restart_times, int global_evap_times)
{

    _past_pheromone = get_pheromone(one_minus_rho, past_trail_restart, past_trail_min, global_restart_times, global_evap_times);
    _local_evap_times = global_evap_times;
}

    public void reinforce(double invert_fitness, double one_minus_rho,
                     double past_trail_restart,
                     double past_trail_min,
                     double trail_max,
                     int global_restart_times, int global_evap_times)
{
    pay_evaporation_debt(one_minus_rho, past_trail_restart, past_trail_min, global_restart_times, global_evap_times);
    _past_pheromone += invert_fitness;
    // _past_pheromone += invert_fitness / _n_child_leaf;
    if (_past_pheromone > trail_max)
        _past_pheromone = trail_max;
}

    public double prob_weight(double alpha, double beta, double one_minus_rho,
                         double past_trail_restart,
                         double past_trail_min,
                         int global_restart_times, int global_evap_times)
{
    double weight;

    weight = prob_weight_without_child_leaf(
        alpha, beta, one_minus_rho,
        past_trail_restart,
        past_trail_min,
        global_restart_times, global_evap_times);
    weight *= _n_child_leaf;

    return weight;
}

    public double prob_weight_without_child_leaf(
    double alpha, double beta, double one_minus_rho,
    double past_trail_restart,
    double past_trail_min,
    int global_restart_times, int global_evap_times)
{
    double pheromone;

    pheromone = get_pheromone(one_minus_rho, past_trail_restart, past_trail_min, global_restart_times, global_evap_times);

    return Math.pow(pheromone, alpha) * Math.pow(_heuristic, beta);
}

    public int choose_child_with_prob(
    double alpha, double beta, double one_minus_rho,
    double past_trail_restart,
    double past_trail_min,
    int global_restart_times,
    int global_evap_times,
    double elite_prob)
{
    // IMPORTANCE NOTE: Remember to check won't visit before go to this function

    double[] weights = new double[2];
    int i;
    double move_prob;

    for (i = 0; i < 2; i++)
        weights[i] = childs.get(i).prob_weight_without_child_leaf(
            alpha, beta, one_minus_rho,
            past_trail_restart, past_trail_min,
            global_restart_times, global_evap_times);

    if (elite_prob != 0.0 && Math.random() < elite_prob)
    {
        if (weights[1] > weights[0])
            return 1;
        else
            return 0;
    }

    for (i = 0; i < 2; i++)
        // weights[i] = childs[i].prob_weight(alpha, beta, one_minus_rho, past_trail_restart, past_trail_min, global_restart_times, global_evap_times);
        weights[i] *= childs.get(i).get_n_child_leaf();

    move_prob = weights[1] / (weights[0] + weights[1]);
    if (Math.random() < move_prob)
        return 1;
    else
        return 0;
}

    public int get_n_child_leaf() { return _n_child_leaf; }
}

class Wont_Visit_Node extends Node_Base<Wont_Visit_Node>{

    public boolean _wont_visit;
    public int _local_wont_visit_restart_times;

void _restart_if_needed(int global_wont_visit_restart_times)
{
    if (_local_wont_visit_restart_times < global_wont_visit_restart_times)
    {
        _wont_visit = false;
        _local_wont_visit_restart_times = global_wont_visit_restart_times;
    }
}

void set_wont_visit(int global_wont_visit_restart_times)
{
    _wont_visit = true;
    _local_wont_visit_restart_times = global_wont_visit_restart_times;
    if (parent != null)
        parent._check_wont_visit(global_wont_visit_restart_times);
}

void _check_wont_visit(int global_wont_visit_restart_times)
{
    if (childs.get(0).get_wont_visit(global_wont_visit_restart_times) && childs.get(1).get_wont_visit(global_wont_visit_restart_times))
        set_wont_visit(global_wont_visit_restart_times);
}

boolean get_wont_visit(int global_wont_visit_restart_times)
{
    _restart_if_needed(global_wont_visit_restart_times);
    return _wont_visit;
}

Wont_Visit_Node(Wont_Visit_Node parent)
{
    this.parent = parent;
    _wont_visit = false;
    _local_wont_visit_restart_times = 0;
}


}

class Building_Node extends Node_Base<Building_Node> {

    double _centroid_x, _centroid_y; // coordinate
    int _n_child_leaf;

double make_heuristic(int city_index)
{
    return compute_heuristic(distance_with_coordinate(city_index, _centroid_x, _centroid_y));
}

int get_n_child_leaf()
{
    return _n_child_leaf;
}


Building_Node(boolean is_leaf)
    : Node_Base<Building_Node>(null, is_leaf) {}

Building_Node(Building_Node *parent, double centroid_x, double centroid_y, int n_child_leaf, boolean is_leaf)
    : Node_Base<Building_Node>(parent, is_leaf)
{
    this._centroid_x = centroid_x;
    this._centroid_y = centroid_y;
    _n_child_leaf = n_child_leaf;
}
}

class Leaf extends public Node {

public int _city_index;

public int get_city_index() { return _city_index; }
public (Node parent, double heuristic, int city_index){
    super(parent, heuristic, 1);
    _city_index = city_index;
        }

Leaf(
    int city_index,
    double heuristic)
    : Node(null, null, false), Leaf_Base(city_index)
{
    _heuristic = heuristic;
}
}

class Building_Leaf : public Building_Node{
Building_Building_Leaf(Building_Node *parent, double centroid_x, double centroid_y, int city_index)
    : Building_Node(parent, centroid_x, centroid_y, 1, true), Leaf_Base(city_index) {}

public int _city_index;

public Leaf_Base(int city_index) { _city_index = city_index; };
public int get_city_index() { return _city_index; };
}