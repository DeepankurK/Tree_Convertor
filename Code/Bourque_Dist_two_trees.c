/* @Copyright  Louxin Zhang, National University of Singapore, 2019 
  * This program computes the Bourque distance and  1 or 2-Bourque Distance 
  * It reads two rooted labeled trees from different files and then calculate
  * the Bourque distance or k-bourque distances and output to standard output
  *  using printf
  *
  *  Tree file:
  *   edges are one by one.
    For instance  "a b_c_d" represents an edge from
  *  a node that is labeled with a to a node that is labeled with 
  *  a subset {b, c, d} of three labels b, c, d; 
      the root is assumed to the the first node.
  *
  *  For example, the following five four lines specify the tree 
      with 4 nodes which are labled with 0, 1, {2, 5}, 3 and 4  and 
      whose root is 0:
  *  0 1
  *  1 3
  *  3 2_5
  *  3 4
  *
  *  LIMITATION: (1) the program uses array as data strucutre to 
      store trees; (2) the max. size of the input is set to be 40
      nodes; (3) If a node is labeled with a set of mut, they are seprared
      with '_'. 

    Compile command: gcc Bourque_Dist_two_trees.c -o Bourque
 *  Run command:   Bourque <tree1_file> <tree2_file> order
 *  where order is an integer 0, 1 or 2.

  *
 */
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#define MAXSIZE  40 

/* deinfe the sise of bipartitie graphs used to compute 1 or 2-BD */
#define INFINITY 199
#define EDGE_NO 1000
#define NODE_NO 70 



 

short  Check_Name(char *node_strings[], short  no_nodes, char *str1){
  short i;

  if (str1==NULL) return -1;
  for (i=0; i<no_nodes; i++) {
	if (strcmp(str1, node_strings[i])==0) return i;
  } 
  return -1;
}




short  min_of_two(short  a, short b){
  if (b==0) return a; 
  else if (a==0)  return b;
  
  if (a< b) return a; else return b;
}



void Adjacent_matrix(short no1, short start1[], short end1[], short adjacent[][MAXSIZE]){
 short i, j, k;

 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) {adjacent[i][j]=0;}
 } 

 for (i=0; i<no1; i++) {
     for (k=0; k<no1-1; k++) {	   
       if (start1[k]==i) { adjacent[i][end1[k]]=1; } 
       if (end1[k]==i) { adjacent[i][start1[k]]=1; }  
     }
 }


} /* Adjacent */


void Distance_Comput(short no1, short start1[], short end1[], short distance[][MAXSIZE]){

 short i, j, k;
 short adjacent[MAXSIZE][MAXSIZE];
 short dist;
 short total;

 Adjacent_matrix(no1, start1,  end1,  adjacent);


 total=0;
 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) { 
	   distance[i][j]=adjacent[i][j];
	  if (distance[i][j]!=0) total=1+total; 
   }
 }


 while (total < (no1*no1 - no1)){
/* printf("total[%d]\n", total); */
  for (i=0; i<no1; i++) {
     for (j=0; j<no1; j++) { 
       if (i!=j && distance[i][j]==0) {
	 for (k=0; k<no1; k++) {
	   if (distance[i][k]!=0 && adjacent[k][j]!=0) 
	       distance[i][j]=min_of_two(distance[i][j], 
			      distance[i][k]+adjacent[k][j]);
	 }
	 if (distance[i][j]!=0) { total=1+total; }
       } 
   } /* j */
  } /* for */
 } /* while */
 /* obtain distance matrix */
}



short Check_Map_Image( short j,  short no1,  short map[]) {
   short i; 
   for (i=1; i<=no1; i++) { if (map[i]==j) return i; }
   return 0;
}




void Comput_map1(short map[], short no1, char *tree1_names[], short no2,  
    char *tree2_names[], short no_common, char *common_names[]){
  short i, j;

  for (i=0; i<no1; i++){
    map[i+1]=-1;
    if (Check_Name(common_names, no_common, tree1_names[i])>=0) {
       for (j=0; j<no2; j++) {
          if (strcmp(tree1_names[i], tree2_names[j])==0) {  map[i+1]=j+1; break;}
       }
    }
  }
} /* end comput_map */




short Get_Index(char *mut_name, short  no_mut, char *mut[]){
    short  k;

    for (k=0; k<no_mut; k++) {
         if (strcmp(mut_name, mut[k])==0) return k;
    }
    return -1;
}

void Decode(char *str, char a, short  *num, char *names[]){
   short len, len1;
   short i, j,k;
   char temp[20];


   k=0;
   len=strlen(str);
   len1=0;
   for (i=0; i<len; i++) {
     if (str[i]!='_') { temp[len1]=str[i]; len1 +=1;}
     if (str[i]=='_' || i==len-1) {
            temp[len1]='\0';
            names[*num]=(char *)malloc(len1); strcpy(names[*num], temp);
            *num +=1;
            len1=0;
     }

   } /* for */
}

/* assume the root is the r_th root */
void   Edge_Partitions_Rooted(short no1,  short start1[], short end1[], 
      char *tree1_names[],  short distance[][MAXSIZE],  short no_mut, 
      char *mut_names[], short partition1[][MAXSIZE], short r){
 short i, j, k;
 short  dist;
 short  total;
 short  head, tail;
 short no_mut1;
 char *node_mut[MAXSIZE];
 short ind;



 for (i=0; i<no1-1; i++){
    partition1[i][0]=no_mut;
      if (distance[r][start1[i]] < distance[r][end1[i]]){
	      head=start1[i]; tail=end1[i]; /* start is more close root. */
      } else { head=end1[i]; tail=start1[i];}


   for (j=0; j<no1; j++) {

            no_mut1=0;
            Decode(tree1_names[j], '_', &no_mut1, node_mut);
      if (distance[head][j]<distance[tail][j]){
	    /* partition1[i][j+1]=0; */
            for (k=0; k<no_mut1; k++) {
               ind=Get_Index(node_mut[k], no_mut, mut_names);
                partition1[i][ind+1]=0;
            }
      } else {
          /* partition1[i][j+1]=1; */  /* all the nodes below an edge */
            for (k=0; k<no_mut1; k++) {
               ind=Get_Index(node_mut[k], no_mut, mut_names);
                partition1[i][ind+1]=1;
            }
      }
   } /* j loop */

 } /* i loop */

} /* edge_partition_rooted */


short BD_Rooted_100(short no1, short no2, short p1[][MAXSIZE], 
      short p2[][MAXSIZE],  short map[]){
short d, i, j, k;
short size1, size2, size0, size3;
short  no_mut1, no_mut2;
short  flag0, flag1;
short  common;
short parti_equal_inform[MAXSIZE][MAXSIZE];
short parti_1_equal_inform[MAXSIZE][MAXSIZE];
short second_term;
short k2, j0, t;
short r; 

   no_mut1=p1[0][0]; no_mut2=p2[0][0];
   common=0;
   for (k=1; k<=no_mut1; k++) if  (map[k]>0 ) common=1+common; 

d=0;
for (i=0; i<no1-1; i++) { /* check every edge */
   size1=0; size0=0;
   for (k=1; k<=no_mut1; k++) { if (p1[i][k]==1 && map[k]>0 ) size1 +=1;} 
    /* the size of nodes below the edge */
   for (k=1; k<=no_mut1; k++) { if (p1[i][k]==0 && map[k]>0 ) size0 +=1;} 
   /* the size of nodes not below the edge */
 
      parti_1_equal_inform[i][i]=1;
   for (j=i; j<no1-1; j++) {
      parti_1_equal_inform[i][j]=0;
      size2=0; 
      for (k=1; k<=no_mut1; k++) {
         if (p1[j][k]==1 && p1[i][k]==0 && map[k]>0 ) size2 +=1;
         if (p1[j][k]==0 && p1[i][k]==1 && map[k]>0 ) size2 +=1;
      }

      if (size2==0) { parti_1_equal_inform[i][j]=1;}
   }/* j loop */

   
   for (j=0; j<no2-1; j++) {
      parti_equal_inform[i][j]=0;
      size2=0; 
      for (k=1; k<=no_mut2; k++) { 
	 if (p2[j][k]==1 && Check_Map_Image(k, no_mut1, map)>0) size2 +=1;
      }
      
      if (size1==size2) {
	flag0=0; flag1=0;
        for (k=1; k<=no_mut1; k++) {
         if (map[k]!=-1){
           if (p1[i][k]==1 && p2[j][map[k]]==1) flag1=1+flag1; 
           if (p1[i][k]==0 && p2[j][map[k]]==0) flag0=1+flag0; 
         }
        }
        if (flag1==size1 && flag0==size0 && size0>0 && size1>0) {
             d=d+2; parti_equal_inform[i][j]=1;
        }
      }
   }/* j loop */
  } /* i loop  */

 if (common==no_mut1) return no1-1+ no2-1-d;
 else {
   second_term=0;
   for (i=0; i<no1-1; i++) {
      r=0;
     for (j=0; j<i; j++) { /* edge i has the same part an edge eariler */ 
      if (parti_1_equal_inform[j][i]==1) r=1;
      if (r==1) {break; } 
     } /* j loop */
 
      if (r!=1) {
         k=0;
         for (j=0; j<no2-1; j++) {
            if (parti_equal_inform[i][j]==1) {k=1+k;} 
         }
         if (k>0) {
           k2=0;
           for (t=i; t<no1-1; t++) {
             if(parti_1_equal_inform[i][t]==1) {k2=1+k2;}}
           if (k<k2) second_term=k+second_term; else second_term=k2+second_term;
         }  
       } /* r */
    } /* for i loop  */
    return no1-1+ no2-1 - second_term;
 } /* else if */
	
} /* BD_rooted_100 */





short  H_BD_Rooted3(short no1, short no2, short p1[][MAXSIZE], 
       short p2[][MAXSIZE], short map[]){
  short d, i, j, k;
  short size1, size2, size0;
  short size3;
  short no_mut1, no_mut2;
  short flag0, flag1;
  short common;
  short parti_equal_inform[MAXSIZE][MAXSIZE];
  short parti_1_equal_inform[MAXSIZE][MAXSIZE];
  short map_1t2[MAXSIZE];
  short map_2t1[MAXSIZE];
  short second_term;
  short k2, j0, t;
  short r; 

   
   no_mut1=p1[0][0]; no_mut2=p2[0][0];

   for (i=1; i<=no_mut2; i++) { map_2t1[i]=-1;}

   common=0;
   for (i=1; i<=no_mut1; i++) {
       map_1t2[i]=-1;
       for (j=1; j<=no_mut2; j++) {
           if (map[p1[0][i]] == p2[0][j]) {  
                map_1t2[i]=j;  map_2t1[j]=i;   common+=1;  break; }
       }
   }
   


   d=0;
   for (i=1; i<=no1-1; i++) { /* check every edge */
    size1=0; size0=0;
    for (k=1; k<=no_mut1; k++) { 
       if (map_1t2[k]>0 ) { if (p1[i][k]==1) size1 +=1; else size0 +=1;} 
    } /* size1: the size of nodes below the edge */
 
    parti_1_equal_inform[i][i]=1;
    for (j=i; j<=no1-1; j++) {
      parti_1_equal_inform[i][j]=0;
      size2=0; 
      for (k=1; k<=no_mut1; k++) {
         if (p1[j][k]==1 && p1[i][k]==0 && map[p1[0][k]]>0 ) size2 +=1;
         if (p1[j][k]==0 && p1[i][k]==1 && map[p1[0][k]]>0 ) size2 +=1;
      }

      if (size2==0) { parti_1_equal_inform[i][j]=1; }
   }/* j loop */


   
   for (j=1; j<=no2-1; j++) {
      parti_equal_inform[i][j]=0;
      size2=0; 
      for (k=1; k<=no_mut2; k++) { if (p2[j][k]==1 && map_2t1[k]>0) size2 +=1; }
      
      if (size1==size2) {
	flag0=0; flag1=0;
        for (k=1; k<=no_mut1; k++) {
         if (map_1t2[k]>0){
           if (p1[i][k]==1 && p2[j][map_1t2[k]]==1) flag1=1+flag1; 
           if (p1[i][k]==0 && p2[j][map_1t2[k]]==0) flag0=1+flag0; 
         }
        }
        if (flag1==size1 && flag0==size0 && size0>0 && size1>0) 
           {d=d+2; parti_equal_inform[i][j]=1;}
      }
   }/* j loop */
 } /* i loop  */

 if (common==no_mut1) { 
     return no1-1+ no2-1-d; 
 } else {
   second_term=0;
   for (i=1; i<=no1-1; i++) {
      r=0;
     for (j=1; j<i; j++) { /* edge i has the same part an edge eariler */ 
      if (parti_1_equal_inform[j][i]==1) r=1;
      if (r==1) {  break; } 
     } /* j loop */
 
      if (r!=1) {
         k=0;
         for (j=1; j<=no2-1; j++) {
            if (parti_equal_inform[i][j]==1) {k=1+k;} 
         }
         if (k>0) {
           k2=0;
           for (t=i; t<=no1-1; t++) {
             if(parti_1_equal_inform[i][t]==1) {k2=1+k2;}
           }
           if (k<=k2) second_term=k+second_term; else second_term=k2+second_term;
         }  
       } /* r */
    } /* for i loop  */
    return no1-1+ no2-1 - second_term;
 } /* else if */
	
}/* HBD3 */


/* part1[i] hodl the parttion in the ith subgragh */
/* part1[i][j] holds the jth partition, j>0 in the ith subgraph */
/* part1[i][0] holds the nodes in the subgraph */
/* roote is the r-th node */
void Sub_Partition_Rooted(short parts1[][MAXSIZE][MAXSIZE], short  no1, 
      short start1[], short end1[], char *tree1_names[], 
      short distances1[][MAXSIZE], short order, short r, short no_mut1, 
      char *mut1[], short nbr_size1[]){
  short i,j, k, t, w;
  short num;
  short ind;
  short head;
  short tail;
  short no_local_muts;
  char *local_mut[MAXSIZE];
  short no_nbs;
  short nbrs[MAXSIZE]; 
  short index;
  short nbr_edges[MAXSIZE];
  


  for (i=0; i<no1; i++) {
    no_nbs=0;  /* no of nodes in the order-neighbors of node i */
    no_local_muts=0; /* no of mutations in the neighbore */
    for (j=0; j<no1; j++) {
	if  (distances1[i][j]<= order 
		&& distances1[i][j]+distances1[r][i]==distances1[r][j]) {
	   nbrs[no_nbs]=j;
	   no_nbs =1+no_nbs;
           Decode(tree1_names[j], '_', &no_local_muts, local_mut);
	}
    }
    nbr_size1[i]=no_nbs;

    num=1;
    for (j=0; j<no_local_muts; j++) {
           ind=Get_Index(local_mut[j], no_mut1, mut1);
           parts1[i][0][num]=ind+1;
           num=1+num;
    }
    parts1[i][0][0]=no_local_muts; /* parts1[i][0][0]: No. of muts in the neighbor*/

    ind=0;
    for (j=0; j<no1-1; j++){  /* edge by edges */
       head=start1[j]; tail=end1[j];  /* head is more close to the root r */

       if (distances1[i][head]<=order &&  distances1[i][tail]<=order 
           && distances1[i][head]+distances1[r][i]+1==distances1[r][tail] ) { 
           /* egde is in the neighor and head is more close to i */
           nbr_edges[ind]=j;
           ind=ind+1;
       }
    }
     
    
           
     for (j=0; j<ind; j++) {
          head=start1[nbr_edges[j]]; tail=end1[nbr_edges[j]];
	  for (k=0; k<no_nbs; k++) {
             no_local_muts=0;
             Decode(tree1_names[nbrs[k]], '_', &no_local_muts, local_mut);
             if (distances1[tail][nbrs[k]]+1== distances1[head][nbrs[k]]){
		    /* parts1[i][ind][k]=1; */
                   for (t=0; t<no_local_muts; t++) {
                     index=Get_Index(local_mut[t], no_mut1, mut1);
                    for (w=1; w<=num-1; w++) {
                      if (parts1[i][0][w]==index+1) { parts1[i][j+1][w]=1; }
                    } 
                   }    
	      }  else {  /* parts1[i][ind][k]=0; */
                   for (t=0; t<no_local_muts; t++) {
                    index=Get_Index(local_mut[t], no_mut1, mut1);
                    for (w=1; w<=num-1; w++) {
                      if (parts1[i][0][w]==index+1) { parts1[i][j+1][w]=0; }
                    } 
                   }    
             }
          }
    } /* j for */
  

  } /* i for */
} /* end of subpatition_rooted */



void  Remove_Space(char h_line[], short len){
      short i;

      for (i=0; i<len-1; i++) h_line[i]=h_line[i+1];
      h_line[len-1]='\0';
}

void dijkstra(short G[NODE_NO][NODE_NO], short n, short startnode,
     short distance[], short pred[]) {

        short cost[NODE_NO][NODE_NO];
        short visited[NODE_NO],count,mindistance,nextnode,i,j;

        for(i=0;i<n;i++) {
          for(j=0;j<n;j++) {
            if (j!=i)   cost[i][j]=G[i][j];
            else if (i==j) { cost[i][j]=INFINITY; }
          }
        }

        /*  initialize pred[],distance[] and visited[] */
        for(i=0;i<n;i++) {
          distance[i]=cost[startnode][i]; pred[i]=startnode; visited[i]=0;
        }
      distance[startnode]=0; visited[startnode]=1; count=1;

      while(count<n-1) {
          mindistance=INFINITY;

          /* nexnode gives the node at minimum distance */
          for(i=0;i<n;i++) {
            if((distance[i]< mindistance ) && !visited[i] ){
                mindistance=distance[i]; nextnode=i;
            }
          }


          /* check if a better path exists through nextnode */
          visited[nextnode]=1;
          for(i=0;i<n;i++) {
             if(!visited[i])
                if(mindistance+cost[nextnode][i]<distance[i]) {
                  distance[i]=mindistance+cost[nextnode][i];
                  pred[i]=nextnode;
             }
          }
              count++;
        } /* end of while */

} /* disji */


/* output graph size  no+2*/
/* covered_l[0] is the no of nodes covered by matching
 *  *  * followed by the indices of covered nodes */
void  Construct_ResidualGraph(short G[][NODE_NO], short BG[][4],
     short  edge_no,  short cost_l[], short cost_r[], short no){

     short  i,j,  covered[NODE_NO],   num;
     short s, t;

     s=2*no; t=s+1;
     for (i=0; i<=t; i++) {
        G[i][i]=0;
        for (j=i+1; j<=t; j++) { G[i][j]=INFINITY;G[j][i]=INFINITY; }
     }


     for (i=0;  i<s; i++) covered[i]=0; /* uncvered  */
     for (i=0;  i<edge_no; i++) {
       if (BG[i][3]==1) { covered[BG[i][0]]=1; covered[BG[i][1]]=1; }
     }

    for (i=0; i<no; i++) { if (covered[i]==0) G[s][i]=0; }
    for (i=no; i<s; i++) { if (covered[i]==0) G[i][t]=0; }

     for (i=0; i<edge_no; i++) {
        if (BG[i][3]==1) /* matching edge */
           G[BG[i][1]][BG[i][0]]=0;
        if (BG[i][3]==0) /* no-matching edges */
           G[BG[i][0]][BG[i][1]]=BG[i][2]+cost_l[BG[i][0]]-cost_r[BG[i][1]];
     }
} /* end contrstruction of residual graph */



/* modified matches from a shortest path found */
void  update_node_costs(short distance[], short cost_l[], short cost_r[], 
      short no){
    short i;

    for (i=0; i<no; i++){
       cost_l[i]=cost_l[i]+distance[i];
       cost_r[i+no]=cost_r[i+no]+distance[no+i];
    }
} /* update_node_costs */


void     update_BG( short pred[], short BG[][4], short edge_no, short no){
    short v1, v2; /* (v1, v2) is an edge */
    short i, s;

      s=2*no; v2=pred[s+1];
      do {
           v1=pred[v2];
        if (v2 >=no) {
           for (i=0; i<edge_no; i++) {
               if (BG[i][0]==v1 && BG[i][1]==v2) {BG[i][3]=1; break;}
           }
        } else {
           if (v1!=s){
             for (i=0; i<edge_no; i++) {
               if (BG[i][0]==v2 && BG[i][1]==v1) {BG[i][3]=0; break;}
             }
           }
        }
        v2=v1;
      } while(v2!=s);
} /* update BG */



/* BG[i] has left node BG[i][0], right node BG[i][1], and cost
  BG[i][2]; BG[i][3]=0 if it is not a matching edge and
  1 otherwise
*/
short Min_Matching(short BG[][4], short no, short edge_no){
    short i, j, k;
    short cost_l[NODE_NO], cost_r[NODE_NO];
    short G[NODE_NO][NODE_NO]; /* max value of dim is 30 */
    short size, weight;
    short distance[NODE_NO], pred[NODE_NO];

  /* edge always from left to right */
   for (i=0; i<no; i++) { cost_l[i]=0; cost_r[no+i]=0;}
   for (i=0; i<edge_no; i++) { BG[i][3]=0; }
   /* no mathcing edges */
   size=0;
   while (size < no) {
     Construct_ResidualGraph(G, BG, edge_no, cost_l, cost_r, no);
     dijkstra(G,2*no+2, 2*no,  distance,  pred);
     update_node_costs(distance, cost_l, cost_r, no);
     update_BG(pred, BG, edge_no, no);
     size +=1;
   }

   weight=0;
   for (i=0; i<edge_no; i++){
      if (BG[i][3]==1)  { weight=weight+BG[i][2]; }
   }
   return weight;
} /*mini_matching  */


short  Compute_root(short no_edges, short start[], short end[]){
    short i, j;
    short *in;
    short *out;
   

    in=(short *)malloc((no_edges+1)*sizeof(short));
    out=(short *)malloc((no_edges+1)*sizeof(short));
 
    for (i=0; i<no_edges+1; i++) { in[i]=0; out[i]=0; }
    for (i=0; i<no_edges; i++) { in[end[i]] +=1; out[start[i]] +=1; }

    for (i=0; i<no_edges+1; i++) { 
      if (in[i]==0 && out[i]>0)  { j=i; break; } 
    }
     free(in); free(out);
    return j;

}

/*  compile command: gcc HBD_rooted_v7.c -o Bourque
 *  Run command:   Bourque tree1_file tree2_file order
 *  where order is an integer 0, 1 or 2.
 */
void main(int argc, char *argv[]){
FILE *In;
char file_name[100];

short i;
short  no1, no_edges1, no2, no_edges2; 
short  u1, u2;
short  start1[MAXSIZE], end1[MAXSIZE];
short  start2[MAXSIZE], end2[MAXSIZE];
char *tree1_names[MAXSIZE], *tree2_names[MAXSIZE];
char *common_names[MAXSIZE];
char str1[20], str2[20];
short no_common;
short  partit1[MAXSIZE][MAXSIZE];
short  partit2[MAXSIZE][MAXSIZE];
short  distances1[MAXSIZE][MAXSIZE], distances2[MAXSIZE][MAXSIZE]; 
short  map[MAXSIZE];
short  order;
short  parts1[MAXSIZE][MAXSIZE][MAXSIZE];
short  parts2[MAXSIZE][MAXSIZE][MAXSIZE];
short nbr_size1[MAXSIZE], nbr_size2[MAXSIZE];
/*
short code[40];
*/
short  no_mut1, no_mut2;
char  *mut1[MAXSIZE], *mut2[MAXSIZE];

/*
short  j, k;
 */
short  BG[EDGE_NO][4];

short  lp1, lp2;
short ind, BD_value;
short  r1, r2;


   if (argc !=4 ) { 
	   printf("Command: ./a.out tree1_file_name tree2_file name order\n"); 
           exit(10);   
  }

   /* pre-processing: read tree 1 from file 1 */
   no1=0; no_edges1=0;

   In=fopen(argv[1], "r");
   if (In ==NULL) printf("Tree_file_name is not readable\n");

   while (fscanf(In, "%s %s\n", str1, str2)!=EOF){
      u1=Check_Name(tree1_names, no1, str1);
	  if (u1==-1) {
		u1=no1;
		tree1_names[no1]=(char *)malloc(strlen(str1)+1);
		strcpy(tree1_names[no1], str1);
		no1=1+no1;
	  } /* the nodes are there */

	  u2=Check_Name(tree1_names, no1, str2);
	  if (u2==-1) {
		u2=no1,  no1=1+no1;
		tree1_names[u2]=(char *)malloc(strlen(str2)+1);
		strcpy(tree1_names[u2], str2);
	  }

      start1[no_edges1]=u1; end1[no_edges1]=u2; no_edges1 +=1;
      /* use two array to store an edges */
   }
   fclose(In);

   no_mut1=0;
   for (i=0; i<no1; i++) {
      Decode(tree1_names[i], '_', &no_mut1, mut1);
   }


   r1=Compute_root(no_edges1, start1, end1);


   In=fopen(argv[2], "r");
   if (In ==NULL) printf("Tree_file_name is not readable\n");

   no2=0; no_edges2=0;
   while (fscanf(In, "%s %s\n", str1, str2)!=EOF){
      u1=Check_Name(tree2_names, no2, str1);
	  if (u1==-1) {
		u1=no2;
		tree2_names[no2]=(char *)malloc(strlen(str1)+1);
		strcpy(tree2_names[no2], str1);
		no2=1+no2;
	  } /* the nodes are there */

	  u2=Check_Name(tree2_names, no2, str2);
	  if (u2==-1) {
		u2=no2,  no2=1+no2;
		tree2_names[u2]=(char *)malloc(strlen(str2)+1);
		strcpy(tree2_names[u2], str2);
	  }

      start2[no_edges2]=u1; end2[no_edges2]=u2; no_edges2 +=1;
      /* use two array to store an edges */
   }
   fclose(In);

   no_mut2=0;
   for (i=0; i<no2; i++) {
      Decode(tree2_names[i], '_', &no_mut2, mut2);
   }


   r2=Compute_root(no_edges2, start2, end2);
 
   no_common=0;
   for (i=0; i<no_mut1; i++) {
         if (Get_Index(mut1[i],  no_mut2, mut2)>=0) { 
	    common_names[no_common]= (char *)malloc(strlen(mut1[i])+1);
		strcpy(common_names[no_common], mut1[i]);
		no_common=1+no_common;
         }
    }


   Distance_Comput(no1, start1, end1, distances1);
   Edge_Partitions_Rooted(no1, start1, end1, tree1_names,  distances1,
    no_mut1, mut1, partit1, r1);

   Distance_Comput(no2, start2, end2, distances2);
   Edge_Partitions_Rooted(no2, start2, end2, tree2_names, distances2, 
   no_mut2, mut2, partit2, r2);

   Comput_map1(map, no_mut1, mut1, no_mut2, mut2, no_common,common_names); 


   order=atoi(argv[3]);
   if (order==0) {
   printf("Bourque_Dist: %d\n", BD_Rooted_100(no1, no2, partit1, partit2,  map)); 
	   exit(100);
   }

  Sub_Partition_Rooted(parts1, no1, start1, end1, tree1_names, distances1, order, r1, no_mut1, mut1, nbr_size1);
  Sub_Partition_Rooted(parts2, no2, start2, end2, tree2_names, distances2, order, r2, no_mut2, mut2, nbr_size2);
 
   if (no1<=no2) {
       ind=0;
       for (lp1=0; lp1<no2; lp1++) {
         for (lp2=0; lp2<no2; lp2++) {
          BG[ind][0]=lp1; BG[ind][1]=no1+lp2;
	  if (lp1<no1) BG[ind][2]=H_BD_Rooted3(nbr_size1[lp1], nbr_size2[lp2], 
              parts1[lp1], parts2[lp2],  map);
	  else { BG[ind][2]=nbr_size2[lp2]-1; 
          }
          ind +=1;
         }
      }
      BD_value=0; BD_value=Min_Matching(BG, no2, no2*no2);
      printf("%d-Bourque distance: %d\n", order, BD_value); 
    } /* if */

   if (no1> no2) {
       ind=0;
       for (lp1=0; lp1<no1; lp1++) {
         for (lp2=0; lp2<no1;  lp2++) {
          BG[ind][0]=lp1; BG[ind][1]=no1+lp2;
	  if (lp2<no2) BG[ind][2]=H_BD_Rooted3(nbr_size1[lp1], nbr_size2[lp2], 
             parts1[lp1], parts2[lp2],  map);
	  else BG[ind][2]=nbr_size1[lp1]-1;
          ind +=1;
         }
      }
      BD_value=0; BD_value=Min_Matching(BG, no1, no1*no1);
      printf("%d-Bourque distance: %d\n", order, BD_value); 
    } /* if */
 
    

}
