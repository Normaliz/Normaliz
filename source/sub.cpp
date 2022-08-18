template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::lift_points_to_this_dim(list<vector<IntegerRet> >& Deg1Proj) {
    if (Deg1Proj.empty())
        return;
    size_t dim = Deg1Proj.front().size() + 1;
    size_t dim1 = dim - 1;

    /*if (dim == EmbDim) {
        InEqus.debug_print('+');
    }*/

    list<vector<IntegerRet> > Deg1Lifted;  // to this dimension if < EmbDim   
    
    size_t max_nr_per_thread = 1000000 / omp_get_max_threads();

    size_t nr_to_lift = Deg1Proj.size();
    NrLP[dim1] += nr_to_lift;
    size_t already_lifted = 0;

    bool not_done = true;
    bool has_poly_equs = (PolyEquations.size() > 0);
    bool has_poly_inequs = (PolyInequalities.size() > 0);

    while (not_done) {
        // cout << "Durchgang dim " << dim << endl;

        not_done = false;
        bool message_printed = false;

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();
        

#pragma omp parallel
        {
            int tn;
            if (omp_get_level() == omp_start_level)
                tn = 0;
            else
                tn = omp_get_ancestor_thread_num(omp_start_level + 1);

            size_t nr_points_in_thread = 0;

            size_t ppos = 0;
            auto p = Deg1Proj.begin();
#pragma omp for schedule(dynamic)
            for (size_t i = 0; i < nr_to_lift; ++i) {
                if (skip_remaining)
                    continue;

                for (; i > ppos; ++ppos, ++p)
                    ;
                for (; i < ppos; --ppos, --p)
                    ;

                if ((*p)[0] == 0)  // point done
                    continue;

                if (!not_done && verbose) {
#pragma omp critical
                    {
                        if (!message_printed)
                            verboseOutput() << "Lifting to dimension " << dim << endl;
                        message_printed = true;
                    }
                }

                not_done = true;

#pragma omp atomic
                already_lifted ++;

                try {
                    IntegerRet MinInterval = 0, MaxInterval = 0;  // the fiber over *p is an interval -- 0 to make gcc happy
                    fiber_interval(MinInterval, MaxInterval, *p);
                    // cout << "Min " << MinInterval << " Max " << MaxInterval << endl;
                    IntegerRet add_nr_Int = 0;
                    if (MaxInterval >= MinInterval)
                        add_nr_Int = 1 + MaxInterval - MinInterval;
                    long long add_nr = convertToLongLong(add_nr_Int);
                    if (dim == EmbDim && count_only && add_nr >= 1 && !primitive
                           && Congs.nr_of_rows() == 0 && Grading.size() == 0 && PolyEquations.size() == 0
                       && PolyInequalities.size() == 0){
#pragma omp atomic
                        TotalNrLP += add_nr;
                    }
                    else {  // lift ppoint
                        for (IntegerRet k = MinInterval; k <= MaxInterval; ++k) {
                            INTERRUPT_COMPUTATION_BY_EXCEPTION

                            vector<IntegerRet> NewPoint(dim);
                            for (size_t j = 0; j < dim1; ++j)
                                NewPoint[j] = (*p)[j];
                            NewPoint[dim1] = k;

                            /*if(primitive){ // in this case we must check equations and true inequalities
                                if(InEqusByDim[EmbDim].nr_of_rows() > 0){
                                    if(!v_non_negative(InEqusByDim[dim].MxV(NewPoint)))
                                        continue;
                                }
                            }*/

                            if(has_poly_equs && !PolyEquations.check(NewPoint, true, true)) // true = equations, true  = exact length
                                continue;
                            if(has_poly_inequs && !PolyInequalities.check(NewPoint, false, true)){ // false = inequalities, true  = exact length
                                continue;
                            }

                            if (dim == EmbDim) {
                                finalize_latt_point(NewPoint, tn); 
                            }
                            else{
                                Deg1Thread[tn].push_back(NewPoint);
                            }
                        }
                    }

                    (*p)[0] = 0;  // mark point as done
                    if (dim < EmbDim)
                        nr_points_in_thread += add_nr;
                    if (nr_points_in_thread > max_nr_per_thread) {  // thread is full
                        skip_remaining = true;
#pragma omp flush(skip_remaining)
                    }

                } catch (const std::exception&) {
                    tmp_exception = std::current_exception();
                    skip_remaining = true;
#pragma omp flush(skip_remaining)
                }

            }  // lifting

        }  // pararllel

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

        for (size_t i = 0; i < Deg1Thread.size(); ++i)
            Deg1Lifted.splice(Deg1Lifted.begin(), Deg1Thread[i]);

        if (dim == EmbDim){
            collect_results(Deg1Lifted);
        }

        // cout << nr_to_lift << " " << already_lifted << endl;
        if(already_lifted == nr_to_lift){
            if(dim1 <= 1){
                if(DoneWithDim.size() > 1)
                    DoneWithDim[1] = true;
                DoneWithDim[0] = true;
            }

            if(dim1 >=1 && DoneWithDim[dim1-1]){
                if(verbose && !DoneWithDim[dim1])
                    verboseOutput() << "Done with dim " << dim1 << " LatticePoints " << NrLP[dim1] << endl;
                DoneWithDim[dim1] = true;
            }
        }
        lift_points_to_this_dim(Deg1Lifted);
        Deg1Lifted.clear();

    }  // not_done

    // cout << "RET FROM " << dim1 << endl;

    if(verbose && dim == EmbDim){
        verboseOutput() << "Complete lattice points so far " << TotalNrLP << endl;
    }

    return;
