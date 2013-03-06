#include "flux.h"
#include "List.h"


// create node 2 node hash table
// input tri[ntri][3] connectivity array
int getn2n(int nn, int ntri, int **tri, List **nhash)
{
  int i,t,k,m,n;
  int n0,n1,n2;


  // node 2 node connectivity
  for (t = 0; t < ntri; t++)
  {
    // loop over 3 nodes of triangle
    for (k = 0; k < 3; k++)
    {
      // get nodes of triangle
      n0 = tri[t][k];
      n1 = tri[t][(k+1)%3];
      n2 = tri[t][(k+2)%3];

      // triangle c connectivity
      //   n2--n1
      //    \ /
      //    n0

      // add self and 2 connected nodes
      nhash[n0]->Ordered_List(n0);
      nhash[n0]->Ordered_List(n1);
      nhash[n0]->Ordered_List(n2);
    }
  }

#if 0
  // Print nhash == n2n
  for (n=0; n < nn;n++)
  {
    printf("\nnode: %2d with %2d nodes, connects to nodes ", n, nhash[n]->max);
    // loop thru nodes connected to node n
    for (i=0; i < nhash[n]->max; i++)
      printf(" %d", nhash[n]->list[i]);
  }
  printf("\n\n");
#endif

  return 0;
}


// create CRS arrays from n2n hash table
int getCRS(int nn, List **nhash, int **IA, int **JA, int **IAU)
{
  int i,m,n,t,k;
  int n0,n1,n2;
  int max;

  // allocate IA[nn+1]
  if (((*IA) = (int*)malloc((nn+1)*sizeof(int))) == NULL){
    printf("\nCould not allocate memory for cdt"); exit(0);
  }
  // allocate IAU[nn]
  if (((*IAU) = (int*)malloc((nn)*sizeof(int))) == NULL){
    printf("\nCould not allocate memory for cdt"); exit(0);
  }


// ========== IA = element number that starts each row ============


  // IA[nn] == elem-# that starts row n
  int nnz = 0;
  for (n=0; n < nn;n++)
  {
    (*IA)[n] = nnz;
    nnz += nhash[n]->max;
  }
  (*IA)[nn] = nnz;


// ========== IAU = element number on diagonal of each row ========


  for (n=0; n < nn;n++)
    (*IAU)[n] = (*IA)[n] + nhash[n]->Index(n);


// ========== JA = all elems of nhash listed in order =============

  // allocate for JA[nnz]
  if (((*JA) = (int*)malloc((nnz)*sizeof(int))) == NULL){
    printf("\nCould not allocate memory for cdt");
    exit(0);
  }

  // Create JA
  k = -1;
  for (n=0; n < nn;n++)
  {
    for (i=0; i < nhash[n]->max; i++)
    {
      m = nhash[n]->list[i];
      (*JA)[++k] = m;
    }
  }


// ==================== Print Zone ================================

#if 0
  // print IA[nn+1]
  for (n=0; n <= nn;n++)
    printf("IA[%2d] = %2d\n",n,(*IA)[n]);
#endif

#if 0
  // Print IAU[nn]
  for (n=0; n < nn;n++)
    printf("IAU[%2d] = %2d\n",n, (*IAU)[n]);
#endif

#if 0
  // Print JA[nnz]
  for (k=0; k < nnz;k++)
      printf("JA[%2d] = %2d\n", k, (*JA)[k]);
#endif

  return (nnz);
}











#if 0
void make_nbrs(int nn, int ntri, int tri[][3], int nbr[][3])
{
  int i,j,k,c;
  int s,m,n;
  int n0,n1,n2;

  // Initialize nbr[ntri][3] to -1
  for (c = 0; c < ntri; c++)
    for (s = 0; s < 3; s++)
      nbr[c][s] = -1;

  List ** nhash;
  nhash = new List*[nn];

  for (n=0; n < nn;n++)
    nhash[n] = new List();

  // node to element connectivity == nhash[nn]->[*]
  for (c = 0; c < ntri; c++)
    for (s = 0; s < 3; s++)
    {
      n = tri[c][s];
      nhash[n]->Add_To_List(c);
    }

  // make neighbors nbr[ntri][3]
  for (c = 0; c < ntri; c++)
  {
    for (s = 0; s < 3; s++)
    {
      switch(s)
      {
        case 0: n0 = tri[c][0]; n1 = tri[c][1]; break;
        case 1: n0 = tri[c][1]; n1 = tri[c][2]; break;
        case 2: n0 = tri[c][2]; n1 = tri[c][0]; break;
      }
      for (i=0; i < nhash[n0]->max && nbr[c][s] < 0; i++)
      {
        m = nhash[n0]->list[i];
        if (m == c)//dont add self as neighbor
          continue;
        if (nhash[n1]->Is_In_List(m))
          nbr[c][s] = m;
      }
    }
  }


  printf("\n====================================\n");
  n = 0;
  printf("checking max elems for node, n = %d\n",n);
  printf("max = %d\n",nhash[n]->max);

  c = 6;
  printf("adding triangle = %d\n", c);
  nhash[n]->Add_To_List(6);

  printf("\nchecking max elems for node, n = %d\n",n);
  printf("max = %d\n",nhash[n]->max);


  printf("\nchecking last value added\n",n);
  c = nhash[n]->list[2];
  printf("last value = %d\n",c);

  printf("\nfinding index of last value added\n",n);
  c = nhash[n]->Index(6);
  printf("index = %d\n",c);


  printf("\nchecking if 6 is in list\n",n);
  c = nhash[n]->Is_In_List(6);
  printf("is_in_list = %d\n",c);

  printf("\nchecking how many times 6 is in list\n",n);
  c = nhash[n]->Times_In_List(2);
  printf("times in list = %d\n",c);

  c = 6;
  printf("\nadding triangle = %d\n", c);
  nhash[n]->Add_To_List(6);

  printf("\nchecking max elems for node, n = %d\n",n);
  printf("max = %d\n",nhash[n]->max);

  printf("\nchecking how many times 6 is in list\n",n);
  c = nhash[n]->Times_In_List(6);
  printf("times in list = %d\n",c);

  c = nhash[n]->Check_List(5);
  printf("\nadding 5 if not in list\n");
  printf("max = %d\n",nhash[n]->max);

  c = nhash[n]->Check_List(7);
  printf("\nadding %d if not in list\n", 7);
  printf("max = %d\n",nhash[n]->max);

  c = nhash[n]->Check_List(7);
  printf("\nadding %d if not in list\n", 7);
  printf("max = %d\n",nhash[n]->max);

  for (i=0; i < nhash[n]->max; i++)
    {
      m = nhash[n]->list[i];
      printf("nhash[n]->list[%d] = %d\n", i, m);
    }

  c = 1;
  while (nhash[n]->Is_In_List(c))
    nhash[n]->Delete_From_List(c);

  printf("\ndeleting %d from list\n", c);
  printf("max = %d\n",nhash[n]->max);

  for (i=0; i < nhash[n]->max; i++)
    {
      m = nhash[n]->list[i];
      printf("nhash[n]->list[%d] = %d\n", i, m);
    }

  printf("\nreplaced all zeros with ones\n");
  nhash[n]->Replace(6,7);

  for (i=0; i < nhash[n]->max; i++)
    {
      m = nhash[n]->list[i];
      printf("nhash[n]->list[%d] = %d\n", i, m);
    }

  FILE * fp = NULL;
  char fname[] = "test.out";
  if ((fp=fopen(fname,"w")) == NULL)
  {
    printf("\nCould not open file <%s>\n",fname);
    exit(0);
  }
  nhash[n]->print(fp);
  //fprintf(fp,"this is a test\n");
  fclose(fp);


  printf("\n====================================\n");
  printf("\n====================================\n");

  for (n=0; n < nn;n++)
    delete nhash[n];
  delete[] nhash;


}
#endif
