
class public class TreeMap {
    static double neighbour_prob, elite_prob;

    vector<Tree_Edge *> _tree_edge_ptrs;
    Wont_Visit_Tree *_wont_visit_tree_ptr;
    int
        _global_wont_visit_restart_times,
        _global_evap_times,
        _global_restart_times,
        _num_city;
    double _past_trail_restart, _past_trail_min;
    

    Tree_Map(const int&num_city, const struct problem &instance) {
        int i;
        Building_Tree building_tree (instance);
        Building_Node * building_root_ptr = building_tree.get_root_ptr();

        this->_num_city = num_city;
        _tree_edge_ptrs.resize(num_city - 1);
        for (i = 0; i < num_city - 1; i++) // Do not go from the end city
        {
            _tree_edge_ptrs[i] = new Tree_Edge(num_city, i, building_root_ptr);
        }

        _wont_visit_tree_ptr = new Wont_Visit_Tree(num_city, building_root_ptr);

        _global_wont_visit_restart_times = 0;
        _global_restart_times = 0;
        _global_evap_times = 0;
    }

    // ~Tree_Map() {
    //     delete _wont_visit_tree_ptr;
    //     for (auto & _ptr :_tree_edge_ptrs)
    //     delete _ptr;
    // }

    void 

    _set_wont_visit(const int&city_index) {
        _wont_visit_tree_ptr -> set_wont_visit(city_index, _num_city, _global_wont_visit_restart_times);
    }

    void 

    reinforce(
    const ant_struct &an_ant,
    double rho,
    double trail_max) {
        int i;
    const double invert_fitness = 1.0 / an_ant.fitness;
    const double one_minus_rho = 1 - rho;

        for (i = 0; i <= an_ant.tour_size - 4; i++) {
            _tree_edge_ptrs[an_ant.tour[i]]->reinforce(
                    an_ant.tour[i + 1],
                    invert_fitness,
                    one_minus_rho,
                    _past_trail_restart,
                    _past_trail_min,
                    trail_max,
                    _global_restart_times,
                    _global_evap_times);
        }
    }

    void 

    _reset_wont_visit() {
        _global_wont_visit_restart_times += 1;
    }

    void 

    global_evaporate(double past_trail_min) {
        _past_trail_min = trail_min;
        _global_evap_times += 1;
    }

    void 

    global_restart(double past_trail_restart) {
        _past_trail_restart = past_trail_restart;
        _global_restart_times += 1;
        _global_evap_times = 0;
    }

    void 

    choose_route(
            ant_struct &an_ant,
    double neighbour_prob,
    double alpha,
    double beta,
    double rho,
            long int&n_tours,
    const int&nn_ants,
            long **nn_list,
    double elite_prob) {
    double one_minus_rho = 1 - rho;
        int current_city, i;

        _reset_wont_visit();
        ant_empty_memory( & an_ant);

        an_ant.tour_size = 1;
        an_ant.tour[0] = 0;
        an_ant.visited[0] = TRUE;
        an_ant.visited[_num_city] = TRUE; // virtual city

        current_city = an_ant.tour[an_ant.tour_size - 1];
        while (current_city != _num_city - 1) {
            current_city = _tree_edge_ptrs[current_city]->choose_next_city(
                    _wont_visit_tree_ptr,
                    neighbour_prob,
                    alpha,
                    beta,
                    one_minus_rho,
                    _past_trail_restart,
                    _past_trail_min,
                    _global_restart_times,
                    _global_evap_times,
                    _global_wont_visit_restart_times,
                    nn_ants,
                    nn_list[current_city],
                    _num_city,
                    elite_prob);
            an_ant.tour[an_ant.tour_size] = current_city;
            an_ant.visited[current_city] = TRUE;
            an_ant.tour_size++;

            _set_wont_visit(current_city);
        }

        an_ant.tour[an_ant.tour_size++] = _num_city;      // virtual city
        an_ant.tour[an_ant.tour_size++] = an_ant.tour[0]; // to form a tour

        for (i = an_ant.tour_size; i < _num_city + 1; i++)
            an_ant.tour[i] = 0;
        an_ant.fitness = compute_fitness(an_ant.tour, an_ant.visited, an_ant.tour_size, an_ant.packing_plan);

        n_tours += 1;
    }

    double 

    leaf_pheromone(
    const int&i,
    const int&j,
    double rho) {
    const double one_minus_rho = 1 - rho;
        return _tree_edge_ptrs[i]->leaf_pheromone(
                j,
                one_minus_rho,
                _past_trail_restart,
                _past_trail_min,
                _global_restart_times,
                _global_evap_times);
    }

    double 

    node_branching(double lambda) {
        int i, j, sum_num_branches = 0;
        double min_pheromone, max_pheromone, cutoff, pheromone;
    const int num_city = instance.n - 1;

        for (i = 0; i < num_city - 1; i++) {
            min_pheromone = 2;
            max_pheromone = -1;
            for (j = 1; j < num_city; j++) {
                if (i == j)
                    continue;
                pheromone = tree_map -> leaf_pheromone(i, j, rho);
                if (pheromone > max_pheromone)
                    max_pheromone = pheromone;
                if (pheromone < min_pheromone)
                    min_pheromone = pheromone;
            }

            cutoff = min_pheromone + lambda * (max_pheromone - min_pheromone);
            for (j = 1; j < num_city; j++) {
                if (tree_map -> leaf_pheromone(i, j, rho) > cutoff)
                    sum_num_branches += 1;
            }
        }

        return (sum_num_branches / (double) (num_city * 2));
    }

    static void tree_map_init() {
        tree_map = new Tree_Map(instance.n - 1, instance);
    }

    static void tree_map_deallocate() {
        delete tree_map;
    }

    static void tree_map_construct_solutions() {
        int i;
        for (i = 0; i < n_ants; i++)
            tree_map -> choose_route(ant[i], neighbour_prob, alpha, beta, rho, n_tours,
                    nn_ants,
                    instance.nn_list, elite_prob);
    }
}