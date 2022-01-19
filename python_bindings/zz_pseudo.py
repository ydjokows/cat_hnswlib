
### Method to change is this one in the 

  .knn_query2( vects, k = topk, conditions = filters   )




    def hnsw_init(dirin="", ndim=512 , efsearch= 32):
        clientdrant = Box({})
        import hnswlib
        p = hnswlib.Index(space   = 'l2', dim = ndim) 
        p.load_index(dirin + "/index.bin")
        p.set_ef( efsearch)   #### Speed/precision        

        df = pd_read_file(dirin + "/index.parquet")
        df = df.set_index('idx').to_dict()
        index  = index['id']   ## 0 --> siid 
        indexg = index.get('genre_id', {})
        clientdrant.p = p ; clientdrant.index = index ; clientdrant.indexg = indexg



    def hnsw_useremb_get_topk(vecti, genreid=[2323,232,3232], topk=100, dimvect=512, filter_cond='must'):
           ####  genrei: list_siid
           global clientdrant
          
           vecti = np.array([ float(x) for x in  vecti.split(",")] ,  dtype='float32')
           idxall = [] 
           vects  = [ vecti ]  #  [ vecti for i in range(0, len(genreid)) ] 
           for gi in genreid :
              filters    = [[ (False, int(gi)) ]]    ### (isexcluded:false, catid) 
              idxallj,_  = clientdrant.p.knn_query( vects, k = topk, conditions = filters   )
              idxall.append( idxallj[0] )

          return idxall 
        
        
    #### New code    
        
    def hnsw_useremb_get_topk(vecti="2330.324,0.34234,0.2324", genreid=[2323,232,3232], topk=100, dimvect=512, filter_cond='must'):
           ####  genrei: list_siid
           global clientdrant
           vecti = np.array([ float(x) for x in  vecti.split(",")] ,  dtype='float32')
           idxall = [] 
           vects  = [ vecti ]  #  [ vecti for i in range(0, len(genreid)) ] 

           idxall  = clientdrant.p.knn_query_new( vects, k = topk, conditions = filters   )
           return idxall 
                
        
        
        
        
        
        
