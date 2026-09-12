// Модульное сито пар (x_i, x_j), i < j: для всех 17 простых (1+x+y) и (1-x-y) — квадраты (или 0) mod p.
// По мотивам six_lines_exact/scripts/sieve.hip.cpp (Codex); изменено: строки i — пачками с хоста, потоки — по j,
// буфер выживших до CAP, бинарный вход/выход. Вход: meta (n k, простые, маски 64 x k), bin (n x k байт).
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#define CHECK(x) do{auto e=(x);if(e!=hipSuccess){fprintf(stderr,"%s\n",hipGetErrorString(e));exit(1);}}while(0)
struct Pair{unsigned int i,j;};
constexpr unsigned long long CAP=64ull*1024*1024;
__global__ void sieve(const unsigned char*res,const unsigned long long*masks,unsigned int n,int k,
                      unsigned int i0,unsigned int i1,Pair*out,unsigned long long*counter){
    unsigned int j=blockIdx.x*blockDim.x+threadIdx.x+i0+1;
    if(j>=n)return;
    unsigned char rj[32];
    for(int l=0;l<k;l++)rj[l]=res[(unsigned long long)j*k+l];
    unsigned int iend=i1<j?i1:j;
    for(unsigned int i=i0;i<iend;i++){
        const unsigned char*ri=res+(unsigned long long)i*k;
        bool ok=true;
        for(int l=0;l<k;l++){ const unsigned long long*m=masks+((size_t)l*128+ri[l])*2; if(!((m[rj[l]>>6]>>(rj[l]&63))&1ULL)){ok=false;break;} }
        if(ok){unsigned long long pos=atomicAdd(counter,1ULL);if(pos<CAP)out[pos]={i,j};}
    }
}
int main(int argc,char**argv){
    if(argc!=4)return 2;
    std::ifstream f(argv[1]);unsigned int n;int k;f>>n>>k;if(!f||k>32)return 3;
    std::vector<int>p(k);for(auto&v:p)f>>v;
    std::vector<unsigned long long>masks(256*k);for(int l=0;l<k;l++)for(int x=0;x<128;x++){f>>masks[(l*128+x)*2];f>>masks[(l*128+x)*2+1];}
    if(!f)return 4;
    std::vector<unsigned char>res((size_t)n*k);
    FILE*fb=fopen(argv[2],"rb");if(!fb)return 5;if(fread(res.data(),1,res.size(),fb)!=res.size())return 6;fclose(fb);
    unsigned char*rr;unsigned long long*mm;Pair*oo;unsigned long long*cc;
    CHECK(hipMalloc(&rr,res.size()));CHECK(hipMalloc(&mm,masks.size()*8));
    CHECK(hipMalloc(&oo,CAP*sizeof(Pair)));CHECK(hipMalloc(&cc,8));
    CHECK(hipMemcpy(rr,res.data(),res.size(),hipMemcpyHostToDevice));
    CHECK(hipMemcpy(mm,masks.data(),masks.size()*8,hipMemcpyHostToDevice));CHECK(hipMemset(cc,0,8));
    hipEvent_t a,b;CHECK(hipEventCreate(&a));CHECK(hipEventCreate(&b));CHECK(hipEventRecord(a));
    const unsigned int CH=256;
    for(unsigned int i0=0;i0+1<n;i0+=CH){
        unsigned int i1=i0+CH<n?i0+CH:n;unsigned int cnt=n-(i0+1);
        hipLaunchKernelGGL(sieve,dim3((cnt+255)/256),dim3(256),0,0,rr,mm,n,k,i0,i1,oo,cc);
    }
    CHECK(hipGetLastError());CHECK(hipEventRecord(b));CHECK(hipEventSynchronize(b));
    float ms;CHECK(hipEventElapsedTime(&ms,a,b));unsigned long long count;CHECK(hipMemcpy(&count,cc,8,hipMemcpyDeviceToHost));
    if(count>CAP){fprintf(stderr,"Overflow %llu\n",count);return 7;}
    std::vector<Pair>out(count);if(count)CHECK(hipMemcpy(out.data(),oo,count*sizeof(Pair),hipMemcpyDeviceToHost));
    FILE*fo=fopen(argv[3],"wb");fwrite(out.data(),sizeof(Pair),count,fo);fclose(fo);
    printf("{\"size\":%u,\"unordered_pairs\":%llu,\"kernel_ms\":%.3f,\"survivors\":%llu}\n",n,(unsigned long long)n*(n-1)/2,ms,count);
    return 0;
}
