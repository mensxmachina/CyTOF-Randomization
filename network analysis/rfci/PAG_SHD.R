pagSHD <- function(PAG1, PAG2){
  
  G1 <- PAG1;
  G2 <- PAG2;
  
  if(is.null(G1) || is.null(G2)){
    return(NA)
  }
  
  nnodes = dim(G1)[1];
  shd = 0;
  for(i in 1:(nnodes-1)){
    
    for(j in (i+1):nnodes){
      
      # o-o
      if(G1[i,j] == 1 && G1[j,i] == 1){
        
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 0;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 1;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 2;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 1;
        }
    
      }
  
      # o->
      if(G1[i,j] == 2 && G1[j,i] == 1){
        
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 1;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 0;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 1;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 2;
        }
        
      }

      # <-o
      if(G1[i,j] == 1 && G1[j,i] == 2){
  
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 1;
        }
        
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 0;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 2;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 2;
        }
        
      }
  
      # <->
      if(G1[i,j] == 2 && G1[j,i] == 2){
      
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 1;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 0;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 1;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 3;
        }
        
      }
      
      # ->
      if(G1[i,j] == 2 && G1[j,i] == 3){
      
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 1;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 0;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 3;
        }
        
      }

      # <-
      if(G1[i,j] == 3 && G1[j,i] == 2){
    
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 1;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 2;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 0;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 3;
        }
        
      }
      
      # 'empty'
      if(G1[i,j] == 0 && G1[j,i] == 0){
        
        # o-o
        if(G2[i,j] == 1 && G2[j,i] == 1){
          shd = shd + 1;
        }
        # o->
        if(G2[i,j] == 2 && G2[j,i] == 1){
          shd = shd + 2;
        }
        # <-o
        if(G2[i,j] == 1 && G2[j,i] == 2){
          shd = shd + 2;
        }
        # <->
        if(G2[i,j] == 2 && G2[j,i] == 2){
          shd = shd + 3;
        }
        # ->
        if(G2[i,j] == 2 && G2[j,i] == 3){
          shd = shd + 3;
        }
        # <-
        if(G2[i,j] == 3 && G2[j,i] == 2){
          shd = shd + 3;
        }
        # 'empty'
        if(G2[i,j] == 0 && G2[j,i] == 0){
          shd = shd + 0;
        }
        
      }
      
    }#for
  }#for
  
  # returning the SHD
  return(shd)
  
}