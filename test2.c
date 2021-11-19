#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define RANGE_MAX 10.0d // максимальное значение для элементов массивов а и с
#define RANGE_MIN -10.0d // минимальное значение для элементов массивов а и с
#define NUM_ELEMENTS 1024 // кол-во эл-в в массивах.
#define THRESHOLD 0.0d
#define EPSILON 1.0e-14d
#ifdef ALIGNED
    #ifdef AVX2
        #define BOUND 32
    #else
        #define BOUND 64 
    #endif
#endif

uint64_t rdtsc()
{
unsigned int lo,hi;
__asm__ __volatile__ ("rdtscp" : "=a" (lo), "=d" (hi)
                               :
                               :"cc", "rcx", "memory");
return ((uint64_t)hi << 32) | lo;
}

int main() {
int i, errorIndex;
double *a, *b1, *b2, *c;  
uint64_t B1_t1, B1_t2;
uint64_t B2_t1, B2_t2;
#ifdef AVX2
double arrayOfTwos[4] __attribute__((aligned(32))) = {2.0d, 2.0d, 2.0d, 2.0d};
double array_TWO_P_FIVE[4] __attribute__((aligned(32))) = {2.5d, 2.5d, 2.5d, 2.5d};
double arrayOf_THRESHOLD[4] __attribute__((aligned(32))) = {THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD};
#else
double arrayOfTwos[8] __attribute__((aligned(64))) = {2.0d, 2.0d, 2.0d, 2.0d, 2.0d, 2.0d, 2.0d, 2.0d};
double array_TWO_P_FIVE[8] __attribute__((aligned(64))) = {2.5d, 2.5d, 2.5d, 2.5d, 2.5d, 2.5d, 2.5d, 2.5d};
double arrayOf_THRESHOLD[8] __attribute__((aligned(64))) = {THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD, THRESHOLD};
#endif

#ifdef ALIGNED
char *tmp_a, *tmp_c, *tmp_b1, *tmp_b2;
uint64_t delta_a, delta_c, delta_b1, delta_b2;
#endif

#ifdef ALIGNED
tmp_a = (char*)malloc(NUM_ELEMENTS * sizeof(double) + (BOUND - 1));
delta_a = (uint64_t)tmp_a - ((uint64_t)tmp_a & (~(BOUND - 1)));
delta_a = (BOUND - delta_a) & (BOUND - 1);

tmp_c = (char*)malloc(NUM_ELEMENTS * sizeof(double) + (BOUND - 1));
delta_c = (uint64_t)tmp_c - ((uint64_t)tmp_c & (~(BOUND - 1)));
delta_c = (BOUND - delta_c) & (BOUND - 1);

tmp_b1 = (char*)malloc(NUM_ELEMENTS * sizeof(double) + (BOUND - 1));
delta_b1 = (uint64_t)tmp_b1 - ((uint64_t)tmp_b1 & (~(BOUND - 1)));
delta_b1 = (BOUND - delta_b1) & (BOUND - 1);

tmp_b2 = (char*)malloc(NUM_ELEMENTS * sizeof(double) + (BOUND - 1));
delta_b2 = (uint64_t)tmp_b2 - ((uint64_t)tmp_b2 & (~(BOUND - 1)));
delta_b2 = (BOUND - delta_b2) & (BOUND - 1);

a = (double*)(tmp_a + delta_a);
c = (double*)(tmp_c + delta_c);
b1 = (double*)(tmp_b1 + delta_b1);
b2 = (double*)(tmp_b2 + delta_b2);
#else
a = (double*)malloc(NUM_ELEMENTS * sizeof(double)); 
c = (double*)malloc(NUM_ELEMENTS * sizeof(double));
b1 = (double*)malloc(NUM_ELEMENTS * sizeof(double));
b2 = (double*)malloc(NUM_ELEMENTS * sizeof(double));
#endif

srand(time(NULL)); 
for (i = 0; i < NUM_ELEMENTS; ++i)
{
    a[i] = RANGE_MIN + (RANGE_MAX - RANGE_MIN) * ((double)rand() / (double)RAND_MAX); 
    c[i] = RANGE_MIN + (RANGE_MAX - RANGE_MIN) * ((double)rand() / (double)RAND_MAX);
    b1[i] = 0.0d;
    b2[i] = 0.0d;
}
B1_t1 = rdtsc();
#ifdef ALIGNED
#pragma vector aligned
#endif
for (i = 0; i < NUM_ELEMENTS; ++i)
{
    if (a[i] > THRESHOLD)
    {
        b1[i] = a[i] + a[i] + c[i];
    } else 
    {
        b1[i] = c[i] * 2.5d + a[i];
    }
}
B1_t2 = rdtsc();

B2_t1 = rdtsc();
#ifdef AVX2
__asm__ __volatile__(".set N, 1024\n\t"
                     "movq $0, %%rcx\n\t"
                     "vmovapd (%[ArrayOfTwos]), %%ymm10\n\t"// load array 2.0
                     "vmovapd (%[Cut_Of]), %%ymm11\n\t"//load array cut_of
                     "vmovapd (%[Two_P_Five]), %%ymm12\n\t"//load array 2.5d
                     "0:\n\t"
                     "vmovupd (%%rdx, %%rcx, 8), %%ymm0\n\t"//load array c[i]
                     "vmovupd (%%rax, %%rcx, 8), %%ymm1\n\t"//load array a[i] 
                     "vmovapd %%ymm1, %%ymm3\n\t"// a[i] -> ymm3 
                     "vcmppd $0xe, %%ymm11, %%ymm1, %%ymm2\n\t"//compare a[i] > cut_of
                     "vfmadd132pd %%ymm10, %%ymm0, %%ymm1\n\t"//fma 2.0d * a[i] + c[i]
                     "vfmadd132pd %%ymm12, %%ymm3, %%ymm0\n\t"//fma 2.5d * c[i] + a[i]
                     "vblendvpd %%ymm2, %%ymm1, %%ymm0, %%ymm1\n\t"
                     "vmovupd %%ymm1, (%%rsi, %%rcx, 8)\n\t" 
                     "addl $4, %%ecx\n\t"
                     "cmpl $N, %%ecx\n\t"
                     "jl 0b\n\t" 
                     :"=c"(i)
                     :"a"(a),"d"(c), "S"(b2), [ArrayOfTwos]"r"(arrayOfTwos), [Cut_Of]"r"(arrayOf_THRESHOLD), [Two_P_Five]"r"(array_TWO_P_FIVE)
                     :"cc","ymm0", "ymm1", "ymm2","ymm3","ymm10","ymm11","ymm12","memory");
#else 
__asm__ __volatile__(".set N, 1024\n\t"
                     "movq $0, %%rcx\n\t"
                     "vmovapd (%[ArrayOfTwos]), %%zmm10\n\t"// load array 2.0
                     "vmovapd (%[Cut_Of]), %%zmm11\n\t"//load array cut_of
                     "vmovapd (%[Two_P_Five]), %%zmm12\n\t"//load array 2.5d
                     "0:\n\t"
                     "vmovupd (%%rdx, %%rcx, 8), %%zmm0\n\t"//load array c[i]
                     "vmovupd (%%rax, %%rcx, 8), %%zmm1\n\t"//load array a[i] 
                     "vmovapd %%zmm1, %%zmm2\n\t"// a[i] -> zmm2 
                     "vcmppd $0xe, %%zmm11, %%zmm1, %%k1\n\t"//compare a[i] > cut_of
                     "knotb %%k1, %%k2\n\t"
                     "vfmadd132pd %%zmm10, %%zmm0, %%zmm1%{%%k1}%{z}\n\t"//fma 2.0d * a[i] + c[i]
                     "vfmadd132pd %%zmm12, %%zmm2, %%zmm0%{%%k2}%{z}\n\t"//fma 2.5d * c[i] + a[i]
                     "vaddpd %%zmm1, %%zmm0, %%zmm0\n\t"
                     "vmovupd %%zmm0, (%%rsi, %%rcx, 8)\n\t"   
                     "addl $8, %%ecx\n\t"
                     "cmpl $N, %%ecx\n\t"
                     "jl 0b\n\t" 
                     :"=c"(i)
                     :"a"(a),"d"(c), "S"(b2), [ArrayOfTwos]"r"(arrayOfTwos), [Cut_Of]"r"(arrayOf_THRESHOLD), [Two_P_Five]"r"(array_TWO_P_FIVE)
                     :"cc","zmm0", "zmm1", "zmm2","zmm10","zmm11","zmm12","k1","k2","memory");
#endif
B2_t2 = rdtsc();

errorIndex = -1;
for (i = 0; i < NUM_ELEMENTS; ++i)
{
    if (fabs(b1[i] - b2[i]) > EPSILON)
    {
        errorIndex = i;
        break;
    }
}

if (errorIndex < 0)
{
    printf("Verification PASSED\n");
} else 
{
    printf("Verification FAILED \nFIRST offending index = %d \nb1[%d] = %f \nb2[%d] = %f \n", 
            errorIndex, errorIndex, b1[errorIndex],errorIndex, b2[errorIndex]);
}	
printf("CPU cycles for all elements in C-version (B1 array): %ld\nCPU cycles for all elements in Assembler-version (B2 array): %ld\n", B1_t2 - B1_t1, B2_t2 - B2_t1);
printf("CPU cycles for element in C-version (B1 array): %f\nCPU cycles for element in Assembler-version (B2 array): %f\n", (double)(B1_t2 - B1_t1)/ NUM_ELEMENTS, (double)(B2_t2 - B2_t1) / NUM_ELEMENTS);
#ifdef ALIGNED
free(tmp_a);
free(tmp_c);
free(tmp_b1);
free(tmp_b2);
#else
free(a);
free(c);
free(b1);
free(b2);
#endif

return 0;
}
