#pragma once
#include "flow.hpp"
#include <fstream>
namespace next5 {
struct Compact {
    Mesh grid;std::vector<I> full,separate;
    Compact(Parent p,int j,int d,int c):grid(p,j,Cuts(5,{d,c,1,true})),full(grid.count),separate(grid.count){}
    void receive(const Mesh& m){
        for(int row=0;row<grid.count;++row){Five x{};
            for(int j=0;j<5;++j)x[j]=grid.axes[j].values[(row/grid.stride[j])%grid.width[j]];
            I value=m.value(x);
            if(m.donor_piece<0)full[row]=std::max(full[row],value);else separate[row]+=value;
        }
    }
    void seal(){
        grid.weights.swap(full);
        for(int j=0;j<grid.count;++j)grid.weights[j]=std::max(grid.weights[j],separate[j]);
        std::vector<I>().swap(separate);
    }
};
inline std::vector<Compact> projection_tables(bool small){
    const auto tree=shapes();std::vector<Compact> tables;
    for(int j=0;j<24;++j)tables.emplace_back(tree[5][j],j,small?4:7,small?3:4);
    return tables;
}
inline void write_projection(const std::string& filename,int sign,I birth,const std::vector<Compact>& tables){
    std::string temporary=filename+".tmp";std::ofstream out(temporary,std::ios::binary);
    auto word=[&](I v){out.write(reinterpret_cast<const char*>(&v),sizeof v);};
    word(0x46544a35000704);word(sign);word(birth);
    for(const auto& t:tables){word(t.grid.count);
        if(t.grid.count)out.write(reinterpret_cast<const char*>(t.grid.weights.data()),sizeof(I)*t.grid.count);}
    out.close();need(out.good(),"projection write failed");need(!std::rename(temporary.c_str(),filename.c_str()),"projection rename failed");
}
inline I read_projection(const std::string& path,int sign,std::vector<Compact>& tables){
    std::ifstream in(path,std::ios::binary);
    auto word=[&](){I v;need(bool(in.read(reinterpret_cast<char*>(&v),sizeof v)),"short projection");return v;};
    need(word()==0x46544a35000704&&word()==sign,"incorrect small-degree projection input");I birth=word();
    for(auto& t:tables){
        need(word()==t.grid.count,"wrong projection shape");t.grid.weights.resize(t.grid.count);
        if(t.grid.count)need(bool(in.read(reinterpret_cast<char*>(t.grid.weights.data()),sizeof(I)*t.grid.count)),"short projection body");
        for(auto z:t.grid.weights)need(0<=z&&z<=5*S,"invalid five-set weight");
        std::vector<I>().swap(t.full);std::vector<I>().swap(t.separate);
    }
    need(in.peek()==EOF,"trailing projection bytes");return birth;
}
} // namespace next5
