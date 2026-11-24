#include "quartet.hpp"
#include <fstream>
int main(int argc,char** argv)try{
    using namespace next5;
    need(argc==4||argc==5,"usage: gradient_probe snapshot requests answers [multiple]");
    mh::Snapshot snapshot(argv[1]);Quartets q(snapshot);q.multiple=argc==5;
    std::ifstream in(argv[2]);std::ofstream out(argv[3]);
    int family,at;
    while(in>>family>>at){
        std::array<Cell,4> x;
        for(auto& t:x)need(bool(in>>t.d>>t.c>>t.high>>t.tail),"truncated slope request");
        out<<q.slope(mh::Family(family),at,x)<<"\n";
    }
    need(in.eof()&&out.good(),"slope probe I/O");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
