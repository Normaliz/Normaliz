#ifndef ETERMS_H
#define ETERMS_H

#include "AnnaUtils.h"

#include <stdint.h>

/*********  (VARIABLE VOID*) FLEXARRAY  *********/

#define GetLen(A)      ((intptr_t)((A)[0]))
//#define SetLen(A,L)    (A)[0] = (void*)(L)
#define SetLen(A,L)    (A)[0] = (L)
#define GetSize(A)     ((intptr_t)((A)[-1]))
//#define SetSize(A,S)   (A)[-1] = (void*)(S)
#define SetSize(A,S)   (A)[-1] = (S)

#define flexarray_malloc(sz,Type) mymalloc((sz+2)*sizeof(void*),Type)+1
#define flexarray_free(p,Type) \
     myfree((GetSize(p)+2)*sizeof(void*),(flexarray)(p-1),Type)
#define flexarray_realloc(p,newsz,Type) \
     myrealloc((flexarray)(p-1), \
     (GetSize(p)+2)*sizeof(void*), (newsz+2)*sizeof(void*), Type)+1

#define PutLast(A,E)       (A)[++(GetLen(A))] = E
#define GetLast(A)         (A)[(GetLen(A))--]
#define MoveLastToNth(A,N) (A)[N] = GetLast(A)

/*********  INT FLEXARRAY  *********/

#define IntsGetLen(I)       (I)[0]
#define IntsGetSize(I)      (I)[-1]
#define IntsSetLen(I,Len)   (I)[0] = Len
#define IntsSetSize(I,Size) (I)[-1] = Size
#define IntsGetLast(I)      (I)[(I[0])--]
#define IntsPutLast(I,B)    (I)[++(I[0])] = B

#define ints_malloc(sz) (ints)mymalloc((sz+2)*sizeof(int),"int")+1
#define ints_free(I)    myfree((IntsGetSize(I)+2)*sizeof(int),(I-1),"int")
#define ints_realloc(I,newsz) (ints)myrealloc(\
   (I-1), (IntsGetSize(I)+2)*sizeof(int), (newsz+2)*sizeof(int), "int")+1

#define IntsMoveLastToNth(I,N)  I[N] = IntsGetLast(I)

ints ints_init(int n);     
ints ints_dup(ints Ints);

/*********  BITSETS  *********/

#define bs_put_n(s,n) \
        do{if (n<=(8*sizeof(shortbitset))) s |= 1<<(n-1);}while (FALSE)
#define bs_get_n(s,n)  (n<=(8*sizeof(shortbitset))) && (s & (1<<(n-1)))
#define bs_del_n(s,n) \
        do{if (n<=(8*sizeof(shortbitset))) s &= ~(1<<(n-1));}while (FALSE)
 
#define bs_init(n) 0
#define bs_free(s) s=0
#define bs_dup(s) s
#define bs_eq(s1,s2) (s1==s2)                           /* TRUE iff s1==s2 */
#define bs_union_and_assign(s1,s2)  (s1 = ((s1)|(s2)))     /* s1 = s1 U s2 */
#define bs_union(s1,s2)  ((s1)|(s2))                    /* returns s1 U s2 */
#define bs_intersection(s1,s2)  ((s1)&(s2))        /* returns s1 inters s2 */
#define bs_contained(s1,s2)   (((s1)&(s2))==(s1))      /* TRUE iff s1 c s2 */
#define bs_disjoint(s1,s2)   (((s1)&(s2))==0)  /* TRUE iff s1, s2 disjoint */

/*********  ETERMS  *********/

#define eterm_size(n) (((2*(n)+4))*sizeof(int))
#define eterm_last(n) 2*(n)+1

#define _ANNA_ETERM_RUM
#ifdef ANNA_ETERM_RUM
typedef struct eterm_rum_aux {
  eterm * slots;
  int eterm_size;
  int top;
} eterm_rum_stack;

extern eterm_rum_stack ETERM_RUM;

void anna_eterm_rum_init(int n);
eterm eterm_rum_malloc(int size);
void eterm_rum_free(eterm p);

#define malloc_eterm(n)   eterm_rum_malloc(ETERM_RUM.eterm_size)+2
#define eterm_free(t)     eterm_rum_free((t)-2)
#else  /*  COCOA_RUM  */
#define malloc_eterm(n)  (eterm)mymalloc(eterm_size(n),"eterm")+2
#define eterm_free(t)     myfree(eterm_size(eterm_get_indetsNo(t)),t-2,"eterm")
#endif

#define eterm_get_indetsNo(t)  (t)[-2]
#define eterm_degree(t)        (t)[-1]

/* TODO: remove */
#define eterm_get_nth(t,n)     (t)[n]

#define SqFr(t)                (t)[0]
#define Indets(t) ((t)+1+eterm_get_indetsNo(t))
#define eterm_get_OccIndNo(t)  (t)[eterm_get_indetsNo(t)+1]

#define SetSqFr(t,sf)             SqFr(t)=sf
#define eterm_set_degree(t,d)     eterm_degree(t)=d
#define eterm_set_indetsNo(t,n)   eterm_get_indetsNo(t)=n;

#define eterm_put_nth(t,n,e)  \
        do {eterm_degree(t)-=t[n]; eterm_degree(t)+=(t[n]=e); \
        if (e) bs_put_n(SqFr(t),n); else bs_del_n(SqFr(t),n);} while (FALSE) 
#define eterm_put0_nth(t,n)  \
        do {eterm_degree(t)-=t[n]; t[n]=0; bs_del_n(SqFr(t),n);} while (FALSE) 
#define eterm_put_non0_nth(t,n,e)  \
        do {eterm_degree(t)-=t[n]; eterm_degree(t)+=(t[n]=e); bs_put_n(SqFr(t),n);} while (FALSE)

#define eterm_free3(t)  eterm_free(t)
/* do not free the term */

eterm eterm_colon(eterm t1, eterm t2);
eterm eterm_init(int n);  /* initialize an eterm */
eterm eterm_dup(eterm t); /* returns a duplicate of t */
eterm eterm_gcd(eterm t1, eterm t2); /* returns gcd(t1,t2) */
eterm eterm_lcm(eterm t1, eterm t2);
void eterm_mult_and_assign(eterm t1, eterm t2);    /* returns t1*=t2 */

coc_bool sp_BigMult(eterm t1, eterm t2);
coc_bool sp_SmallMult(eterm t1, eterm t2);
coc_bool sp_Mult(eterm t1, eterm t2);
coc_bool sp_eq(eterm t1, eterm t2);
coc_bool sp_eterm_coprime (eterm t1, eterm t2);
int eterm_OccIndets(eterm t);
eterm eterm_colon(eterm t1, eterm t2); /* t2 sparse, returns t1 modified */
void eterm_union_and_assign(eterm t1, eterm t2);

#if 0
coc_bool  eterm_divides(eterm t1, eterm t2);
#endif

#define eterm_divides(t1,t2) \
 ( (bs_contained(SqFr(t1),SqFr(t2))) ? (sp_Mult(t2,t1)) : FALSE )

#define eterm_eq(t1,t2) \
 ( (bs_eq(SqFr(t1),SqFr(t2))) ? (sp_eq(t2,t1)) : FALSE )

#define eterm_which_divides(t1,t2) \
     ( eterm_eq(t1,t2) ? ANSWER_BOTH : \
       ( eterm_divides(t1,t2) ? ANSWER_FIRST : \
	 ( eterm_divides(t2,t1) ? ANSWER_SECOND : ANSWER_NC ) ) )
#if 0
#define eterm_divides(t1,t2) \
   ( ( OccIndNo(t1) <= OccIndNo(t2) )  ?  ( sp_Mult(t2, t1) )  :  FALSE )
#define eterm_divides(t1,t2)  (sp_Mult(t2, t1)) 
#endif

#define eterm_coprime(t1,t2) \
 ( (bs_disjoint(SqFr(t1),SqFr(t2)))  ?  (sp_eterm_coprime(t1, t2))  :  FALSE )


/******************  bgrobner & toric  ******************/

eterm eterm_colon_dup (eterm t1, eterm t2);
void eterm_complete(eterm t);
eterm special_eterm_colon_dup (eterm t1, eterm t2);
eterm ivec_pos2eterm (ivec v);
eterm ivec_neg2eterm (ivec v);

#endif  /* ETERMS_H */
