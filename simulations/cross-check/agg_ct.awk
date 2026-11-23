# aggregate friend_ct output: per checkpoint t, sums of n, C, C^2, S, S^2 and run count
!/^#/ && NF==9 { t=$2; R[t]++; sn[t]+=$3; sc[t]+=$4; scc[t]+=$4*$4; ss[t]+=$5; sss[t]+=$5*$5 }
END { for (t in R) printf "t %s runs %d sum_n %.17g sum_C %.17g sum_C2 %.17g sum_S %.17g sum_S2 %.17g\n", t, R[t], sn[t], sc[t], scc[t], ss[t], sss[t] }
