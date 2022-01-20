
#################### Goal
    move the python Loop  INSIDE the C++ code binding.cpp  code
   
    ---> Loop over the condition with only  1 vector      
             (before only 1 condition and loop over the input vectors)

  

  
#########  
binding.cpp  : add new method:    knn_query_new
    
    
  PYBIND11_PLUGIN(hnswlib) {
        py::module m("hnswlib");
  
          .def("knn_query_new", &Index<float>::knnQuery_return_numpy_new, py::arg("data"), py::arg("k")=1, py::arg("num_threads")=-1, 
                                         py::arg("conditions")=std::vector< std::vector<std::vector< hnswlib::tagtype >>  >())
                                         // Triple Vector  [ [[  tagtype  ]]  ]
    
.....
    
    
       py::object knnQuery_return_numpy(py::object input, size_t k = 1, int num_threads = -1, hnswlib::condition_t &conditions = {}) {

        py::array_t < dist_t, py::array::c_style | py::array::forcecast > items(input);
        auto buffer = items.request();

    
 
         
....

    #### Always (buffer.ndim == 1)     
    rows = 1;
    features = buffer.shape[0];
         
         
     // Loop over the condition with only  1 input vector      (before only 1 condition and loop over the input vectors)
    if(normalize==false) {
        ParallelFor(0, conditions, 
                    num_threads, 
                    [&](size_t ncondition, size_t threadId) {
                          
                           //  hnswlib::SearchCondition search_condition = hnswlib::SearchCondition(conditions);
                      
                           std::priority_queue<std::pair<dist_t, hnswlib::labeltype >> 
                           result = appr_alg->searchKnn((void *) items.data(row), k,  (void *) conditions);                      
                           while(!result.empty()){    
    
    
    

                               
                               
#####                               
 #### Binding usage in Python ##############################################################
  
  ##############################################################################################
  ##### Current Python version ################################################################
    def hnsw_init(dirin="", ndim=512 , efsearch= 32):
        clientdrant = Box({})
                
        import hnswlib
        p = hnswlib.Index(space   = 'l2', dim = ndim) 
        p.load_index(dirin + "/index.bin")
        p.set_ef( efsearch)   #### Speed/precision        

        clientdrant.p = p 



    def hnsw_useremb_get_topk(vecti, genreid=[2323,232,3232], topk=100, dimvect=512, filter_cond='must'):
           ####  genrei: list_siid
           global clientdrant
          
           vecti = np.array([ float(x) for x in  vecti.split(",")] ,  dtype='float32')
           idxall = [] 
                               
           ##### Loop is in python                    
           for gi in genreid :
              filters    = [[ (False, int(gi)) ]]    ### (isexcluded:false, catid) 
              idxallj,_  = clientdrant.p.knn_query( vecti, k = topk, conditions = filters   )
              idxall.append( idxallj[0] )

          return idxall 
        
        

                             
        
#################################################################################### ##########################################                 
    #### New Python version      
        
    def hnsw_useremb_get_topk(vecti="2330.324,0.34234,0.2324", genreid=[2323,232,3232], topk=100, dimvect=512, filter_cond='must'):
           ####  genrei: list_siid
           global clientdrant
           vecti = np.array([ float(x) for x in  vecti.split(",")] ,  dtype='float32')
           idxall = [] 

           #### Loop is in C++                  
           idxall  = clientdrant.p.knn_query_new( vecti, k = topk, conditions = filters   )
           return idxall 
                
        
        
        
        
        
        
