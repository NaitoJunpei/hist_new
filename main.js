// size parameters of figures
var x_base = 20;
var width_graph = 780;
var height_spike = 10;
var height_graph = 60;
var height_hist =54;
var max_repeat = 500;
var max_count = 10; // the number of initial binning positions taken for averaging
var res_graph = 200;

var spike_num;
var onset, offset;
var lv, np;

function LoadStart() {
  return new Promise(function(Main){
 document.getElementById("loading").style.visibility="visible";
 Main();
  });
}

function LoadEnd(){
 document.getElementById("loading").style.visibility="hidden";
}

// data initialization
function ResetData() {
  //document.data.spikes.value = "0.049 0.141 0.225 0.274 0.303 0.320 0.336 0.437 0.478 0.496 0.538 0.553 0.562 0.580 0.632 0.633 0.645 0.659 0.663 0.673 0.678 0.700 0.721 0.728 0.750 0.765 0.771 0.792 0.815 0.838 0.853 0.867 0.905 0.923 0.936 0.947 0.990 1.003 1.021 1.052 1.073 1.106 1.112 1.141 1.153 1.170 1.185 1.213 1.215 1.279 1.285 1.338 1.417 1.462 1.587 1.591 1.764 1.871 1.888 1.944 1.965 2.006 2.013 2.021 2.046 2.063 2.105 2.126 2.136 2.164 2.216 2.288 2.291 2.308 2.393 2.424 2.436 2.471 2.559 2.580 2.602 2.643 2.674 2.718 2.762 2.909 2.947 2.979 3.010 3.032 3.039 3.052 3.103 3.157 3.185 3.217 3.242 3.273 3.326 3.337 3.348 3.375 3.392 3.415 3.426 3.449 3.475 3.489 3.545 3.636 3.656 3.669 3.696 3.718 3.749 3.780 3.849 3.862 3.869 3.935 4.084 4.202 4.222 4.249 4.287 4.348 4.365 4.386 4.416 4.428 4.431 4.442 4.459 4.492 4.501 4.520 4.534 4.550 4.595 4.603 4.634 4.639 4.655 4.662 4.678 4.691 4.710 4.725 4.732 4.753 4.766 4.805 4.820 4.868 4.887 4.891 4.931 4.965 4.991 5.012 5.072 5.083 5.111 5.181 5.220 5.277 5.327 5.431 5.494 5.565 5.838 5.863 5.894 6.014 6.086 6.103 6.119 6.137 6.149 6.168 6.186 6.214 6.264 6.278 6.306 6.353 6.414 6.422 6.450 6.517 6.532 6.598 6.666 6.693 6.711 6.743 6.788 6.803 6.838 6.846 6.863 6.876 6.891 6.909 6.952 6.959 6.976 6.996 7.015 7.028 7.039 7.052 7.057 7.092 7.130 7.148 7.165 7.195 7.226 7.230 7.241 7.247 7.275 7.287 7.302 7.311 7.317 7.326 7.340 7.354 7.381 7.407 7.440 7.466 7.517 7.519 7.583 7.645 7.658 7.676 7.689 7.778 7.788 7.832 7.864 7.884 7.973 8.042 8.167 8.523 8.592 8.644 8.724 8.776 8.809 8.842 8.863 8.892 8.965 8.969 8.981 9.001 9.012 9.025 9.055 9.060 9.085 9.099 9.136 9.168 9.206 9.212 9.237 9.256 9.294 9.301 9.309 9.330 9.367 9.397 9.448 9.565 9.630 9.683 9.736 9.783 9.845 9.867 9.914 9.950 9.985";
 document.data.spikes.value = "1.304 1.317 1.455 1.547 1.565 1.603 1.605 1.628 1.665 1.679 1.684 1.743 1.765 1.767 1.773 1.774 1.806 1.832 1.847 1.863 1.878 1.882 1.909 1.923 1.926 1.939 1.972 1.998 2.043 2.046 2.065 2.088 2.094 2.132 2.142 2.177 2.184 2.193 2.215 2.267 2.291 2.307 2.338 2.397 2.433 2.473 2.518 2.537 2.543 2.580 2.581 2.739 2.766 2.799 2.964 3.082 3.368 3.411 3.512 3.582 3.598 3.710 3.875 3.917 4.146 4.231 4.525 4.872 5.004 5.067 5.091 5.201 5.235 5.310 5.417 5.514 5.554 5.589 5.649 5.668 5.764 5.780 5.794 5.829 5.873 5.900 5.907 5.952 5.979 6.035 6.053 6.092 6.141 6.161 6.189 6.252 6.265 6.292 6.336 6.385 6.448 6.491 6.561 6.656 6.790 6.832 6.970 7.017 7.130 7.342 7.370 7.428 7.448 7.464 7.513 7.528 7.632 7.670 7.683 7.705 7.711 7.718 7.768 7.809 7.815 7.824 7.872 7.881 7.918 7.949 7.953 7.979 7.983 8.061 8.138 8.197 8.252 8.271 8.314 8.323 8.522 8.528 8.540 8.569 8.573 8.584 8.628 8.630 8.658 8.711 8.781 8.854 8.865 9.050 9.154 9.695 9.731 9.833 9.889 9.980";
 return 0;
}

// generating data randomly
var MT = new MersenneTwister();
var Alpha = 2.0*Math.PI*MT.next();
var Beta  = 2.0*Math.PI*MT.next();
var Theta = 2.0*Math.PI*MT.next();
var Amp = 0.3+1.2*MT.next();

var SpikeData = new Array(3);

function RandomData() {
    var t1,t2;
    t1=Solve(0.0,Gamma(1.0));
    document.data.spikes.value = Number(t1.toFixed(3));
    var j=1;

    Alpha = 2.0*Math.PI*MT.next();
    Beta = 2.0*Math.PI*MT.next();
    Theta = 2.0*Math.PI*MT.next();
    Amp = 0.3+1.2*MT.next();
    
    var kappa = Math.random() * 1.25 + 0.75;
    while(1){
        t2=t1+Solve(t1,Gamma(kappa));
        if(t2>TIME) break;
        document.data.spikes.value += " " + t2.toFixed(3);
        t1=t2;
        j++;
    }
    return 0;
}
    
var Base=30.0;
var Amplitude=10.0;
var TIME=50.0/3;
var Period=[2.0/Math.PI,1.41421356/Math.PI,0.8989898/Math.PI];

function Rate_integral(prev_time,new_time){
    return Base*(new_time-prev_time) - Amplitude*Period[0]*Amp*( Math.cos(Alpha+new_time/Period[0]/Amp) - Math.cos(Alpha+prev_time/Period[0]/Amp) ) - Amplitude*Period[1]*Amp*( Math.cos(Beta+new_time/Period[1]/Amp) - Math.cos(Beta+prev_time/Period[1]/Amp) ) - Amplitude*Period[2]*Amp*( Math.cos(Theta+new_time/Period[2]/Amp) - Math.cos(Theta+prev_time/Period[2]/Amp) );
}

function Solve(prev_time,interval){
 var boundary = new Array(2);
 var new_interval;
 boundary[0]=0; boundary[1]=0.5/Base;
 while( Rate_integral(prev_time,prev_time+boundary[1]) < interval ){  boundary[1]+=0.5/Base;  }
 
 while( boundary[1]-boundary[0] > Math.pow(10.0,-6.0) ){
  new_interval=0.5*(boundary[0]+boundary[1]);
  if( Rate_integral(prev_time,prev_time+new_interval) > interval ) boundary[1]=new_interval;
  else                boundary[0]=new_interval; 
 }
 new_interval=0.5*(boundary[0]+boundary[1]);
 if(new_interval<Math.pow(10.0,-8.0)) new_interval=Math.pow(10.0,-8.0);
 return new_interval;
}
function Gamma( kappa ){
    var int_kappa=Math.floor(kappa);
    var frac_kappa=kappa-Math.floor(kappa);
 var x_frac,x_int;
 /*integer part*/
 x_int=0;
 for(var i=0;i<int_kappa;i++){
  x_int+=-Math.log(MT.next());
 }
    /*fractional part*/
 if( frac_kappa < 0.01 ) x_frac=0;
 else{
  var b=(Math.exp(1.0)+frac_kappa)/Math.exp(1.0);
  while(1){
   var u=MT.next();
   var p=b*u;
   var uu=MT.next();
   if(p<=1.0){
    x_frac=Math.pow(p,1.0/frac_kappa);
    if(uu<=Math.exp(-x_frac)) break;
   }
   if(p>1.0){
    x_frac=-Math.log((b-p)/frac_kappa);
    if(uu<=Math.pow(x_frac,frac_kappa-1.0)) break;
   }
  }
 }
 return (x_int+x_frac)/kappa;
}

var time_old = new Array();
// produces a Data object, storing the local time . 
var date_obj = new Date();

// main function
function Main() {
  var spike_time = new Array();
  PostData(spike_time);

  spike_num = spike_time.length;
  // planning to implement a sort function
  onset = spike_time[0] - 0.001 * (spike_time[spike_num - 1] - spike_time[0]);
  offset = spike_time[spike_num - 1] + 0.001 * (spike_time[spike_num - 1] - spike_time[0]);
  
  // transform from absolute times to the passage of times
  time_old[0] = new Date().getTime();
  SpikeRaster(spike_time);
  time_old[1] = new Date().getTime();
  DrawGraph_SSOS(spike_time);   // old method new method
  time_old[3] = new Date().getTime();
  DrawGraph_Kernel12(spike_time); // kernel smoother (with and without reflection boundary)
  time_old[5] = new Date().getTime();
  //DrawGraph_BayesNP(spike_time); // Bayes method for non-poisson spike train
  //time_old[5] = new Date().getTime();
  DrawGraph_Bayes(spike_time); // Bayes method
  time_old[6] = new Date().getTime();
  DrawGraph_HMM(spike_time);  // Hidden Markov Model
  time_old[7] = new Date().getTime();
  
  //document.getElementById("time").innerHTML = "<font size='2pt' face='Arial'>Spike Raster : " + (time_old[1]-time_old[0]) + " ms<br>(A) : " + (time_old[2]-time_old[1]) + " ms<br>(B)-(A) : " + (time_old[3]-time_old[2]) + " ms<br>(C) : " + (time_old[4]-time_old[3]) + " ms<br>(D)-(C) : " + (time_old[5]-time_old[4]) + " ms<br>(E) : " + (time_old[6]-time_old[5]) + " ms<br>(F) : " + (time_old[7]-time_old[6]) + " ms</font>";
  document.getElementById("time").innerHTML = "<font size='2pt' face='Arial'>Computation times of the Javascript codes: (A) : " + (time_old[2]-time_old[1]) + " ms; (B)-(A) : " + (time_old[3]-time_old[2]) + " ms; (C) : " + (time_old[4]-time_old[3]) + " ms; (D)-(C) : " + (time_old[5]-time_old[4]) + " ms; (E) : " + (time_old[6]-time_old[5]) + " ms; (F) : " + (time_old[7]-time_old[6]) + " ms</font>";
  
  //DrawGraph(spike_time, SS(spike_time), "SS");  // old method
  //DrawGraph(spike_time, OS(spike_time), "OS");  // new method
  //DrawGraph(spike_time, Kernel(spike_time), "Kernel"); // kernel smoother
  //DrawGraph(spike_time, Kernel(spike_time), "Kernel2"); // kernel smoother with reflection boundary
  //DrawGraph(spike_time, 0, "HMM"); // Hidden Markov Model
}

// processing input data
function PostData(spike_time) {
  var data_text = document.data.spikes.value.replace(/\r?\n/g," ").replace(/^\s+|\s+$/g,"");
  var data_seq = data_text.split(/[^0-9\.]+/);
  //document.data.spikes.value = ""
  for (var i = 0; i < data_seq.length; i++) {
    spike_time[i] = Number(data_seq[i]);
    // document.data.spikes.value += Math.round(spike_time[i]*1000) + " ";
  }
  // sorting
  spike_time.sort(function(a,b){
      if( a < b ) return -1;
      if( a > b ) return 1;
      return 0;
  });
}

function SpikeRaster(spike_time){
 var names = ['SS','OS','Kernel','Kernel2','Bayes','HMM'];
 names.forEach(function(name){
  var wrap = d3.select('#raster_' + name);
  wrap.select("svg").remove();
  var svg = wrap.append("svg").attr("width",800).attr("height",15);
  
  var line = svg.append("line")
        .attr("x1",x_base)
        .attr("y1",height_spike)
        .attr("x2",x_base+width_graph)
        .attr("y2",height_spike)
        .attr("stroke","black")
        .attr("stroke-width",1);
  for (var i = 0; i < spike_num; i++) {
       x = (spike_time[i] - onset) / (offset - onset);
       var line = svg.append("line")
          .attr("x1",x_base+width_graph * x)
          .attr("y1",0)
          .attr("x2",x_base+width_graph * x)
          .attr("y2",height_spike)
          .attr("stroke","black")
          .attr("stroke-width",1);
  }
 });
}

// Shimazaki-Shinomoto, Omi-Shinomoto
function SSOS(spike_time) {
 var binsize;
 var count = new Array();
 var cost_SS, cost_OS, cost_SS_min, cost_OS_min;
 var w_av, av, va;
 var fano;
 var opt_binsize = new Array(); // [0]:SS, [1]:OS
 lv = 0;
 // compute the local variation Lv measuring the irregularity
 for (var i = 0; i < spike_num - 2; i++) {
  var interval = new Array(2);
     interval[0] = spike_time[i + 1] - spike_time[i];
     interval[1] = spike_time[i + 2] - spike_time[i + 1];
     if ((interval[0] + interval[1]) != 0) lv += 3 * Math.pow(interval[0] - interval[1], 2.0) / Math.pow(interval[0] + interval[1], 2.0) / (spike_num - 2);
     else lv += 3.0 / (spike_num - 2);
 }
 if (lv < 1) np = "regular";
 else np = "bursty";
 // vary the number of bins (max 500) 

 var TT = spike_time.concat(spike_time.map(function(element) {
  return element + (offset - onset);
 }));
 for (var bin_num = 1; bin_num < max_repeat; bin_num++) {
  binsize = (offset - onset) / bin_num;
  cost_SS = 0;
  cost_OS = 0;
  for (var cost_count = 0; cost_count < max_count; cost_count++) {
   start = onset + cost_count * (binsize) / max_count;
   end = offset + cost_count * (binsize) / max_count;
   // initialization of the spike count
   for (i = 0; i < bin_num; i++) {
    count[i] = 0;
   }
   //count the number of spikes
   for (i = 0; TT[i] < end; i++) {
    if (TT[i] >= start) {
     count[Math.floor((TT[i] - start) / binsize)]++;
    }
   }
   // computing the mean and variance of the numbers of spikes in a bin
   av = 0;
   va = 0;
   w_av = 0;
   for (i = 0; i < bin_num; i++) {
    if (count[i] > 2) {
     fano = 2.0 * lv / (3.0 - lv);
    } else {
     fano = 1.0;
    }
    w_av += fano * count[i] / bin_num;
    av += count[i] / bin_num;
    va += count[i] * count[i] / bin_num;
   }
   // computing the cost function
   cost_SS += (2.0 * av - (va - av * av)) / (binsize * binsize);
   cost_OS += (2.0 * w_av - (va - av * av)) / (binsize*binsize);
  }
  cost_SS /= max_count;
  cost_OS /= max_count;
  // updates if the cost is smaller
  if (cost_SS < cost_SS_min || bin_num == 1) {
   cost_SS_min = cost_SS;
   opt_binsize[0] = binsize;
  }
  if (cost_OS < cost_OS_min || bin_num == 1) {
   cost_OS_min = cost_OS;
   opt_binsize[1] = binsize;
  }
 }
 return opt_binsize;
}

// kernel smoother
function Kernel(spike_time){
 var width = new Array(50);
 var cost = new Array(width.length);
 var cost_min,min_index;
 for (var i=0; i<width.length; i++) {
  width[i] = (offset - onset) / (i+1);
  cost[i] = KernelCost(spike_time, width[i]);
  if(cost[i]<cost_min || i==0){
   cost_min = cost[i];
   min_index = i;
  }
 }
 return width[min_index];
}

// cost function of the kernel smoother
function KernelCost(spike_time, width) {
 var A = 0; 
 for (var i=0; i<spike_time.length; i++) {
  for (var j=i+1; j<spike_time.length; j++) {
   var x = spike_time[i]-spike_time[j];
   if (x < 5*width) {
    A = A + 2*Math.exp(-x*x/4/width/width) - 4*Math.sqrt(2)*Math.exp(-x*x/2/width/width);
   }
  }
 }
 return (spike_time.length/width + A/width) / 2 / Math.sqrt(Math.PI);
}

function Bayes(spike_time){
 var n=0; // 10^n < x < 10^(n+1)
 if(offset-onset>1){
  while((spike_time[spike_time.length-1]-spike_time[0])>Math.pow(10,n+1)){
   n += 1;
  }
 }else{
  while((spike_time[spike_time.length-1]-spike_time[0])<Math.pow(10,n)){
   n -= 1;
  }
 }
}

function DrawGraph_SSOS(spike_time){
 //SS
 var wrap = d3.select('#graph_SS');
 wrap.select("svg").remove(); // initialization
 var svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);
 
 var opt = new Array();
 opt = SSOS(spike_time);
 var opt_rate_SS = new Array();
 var rate_max = EstimateRate(spike_time, opt[0], opt_rate_SS);
 var x,y,xx,yy;
 for (var i = 0; i < opt_rate_SS.length; i++) {
     x = i * opt[0] / (offset - onset);
     y = opt_rate_SS[i] / rate_max;
     xx = x_base + width_graph * x;
     yy = height_hist * y;
     if (onset + (i + 1) * opt[0] < offset){
      svg.append("rect").attr("x", xx).attr("y", height_graph-yy).attr("width", width_graph * opt[0] / (offset - onset)).attr("height", yy).attr("fill","#87CEFA").attr("stroke","#67AEDA");
     }else{
      svg.append("rect").attr("x", xx).attr("y", height_graph-height_hist * y).attr("width", width_graph - width_graph * x).attr("height", height_hist * y).attr("fill","#87CEFA").attr("stroke","#67AEDA");
     }

 }
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");
 document.getElementById("optimal_SS").innerHTML = "Optimal bin size = <font color=\"red\">" + opt[0].toFixed(2) + "</font>";

 time_old[2] = new Date().getTime();
 
 //OS
 wrap = d3.select('#graph_OS');
 wrap.select("svg").remove(); // initialization
 svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);
 
 var opt_rate_OS = new Array();
 rate_max = EstimateRate(spike_time, opt[1], opt_rate_OS);
 for (var i = 0; i < opt_rate_OS.length; i++) {
     x = i * opt[1] / (offset - onset);
     y = opt_rate_OS[i] / rate_max;
     xx = x_base + width_graph * x;
     yy = height_hist * y;
     if (onset + (i + 1) * opt[1] < offset){
      svg.append("rect").attr("x", xx).attr("y", height_graph-yy).attr("width", width_graph * opt[1] / (offset - onset)).attr("height", yy).attr("fill","#7FFFD4").attr("stroke","#5FDFB4");
     }else{
      svg.append("rect").attr("x", xx).attr("y", height_graph-height_hist * y).attr("width", width_graph - width_graph * x).attr("height", height_hist * y).attr("fill","#7FFFD4").attr("stroke","#5FDFB4");
     }
 }
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");

 document.getElementById("optimal_OS").innerHTML = "Optimal bin size = <font color=\"red\">" + opt[1].toFixed(2) + "</font>&nbsp;&nbsp;&nbsp;&nbsp;Irregularity is estimated as Lv = <font color=\"red\">" + lv.toFixed(2) + "</font>";
}

function DrawGraph_Kernel12(spike_time){
 // Kernel(C)
 var wrap = d3.select('#graph_Kernel');
 wrap.select("svg").remove(); // initialization
 var svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);
 
 //var opt = Kernel(spike_time);
 var opty1 = new Array();
 var opty2 = new Array();
 //var maxy = kern12(spike_time, opt, opty1, opty2);
 var res = kernel_rate(spike_time,opty1,opty2);
 var maxy = res[0];
 var opt = res[1];
 var xy1 = new Array();
 for (var i = 0;i<opty1.length;i++) {
  xy1[i] = [x_base + Math.round(i*width_graph/(opty1.length-1)), height_graph - Math.round(height_graph*opty1[i]/(1.2*maxy))];
 }
 xy1.unshift([x_base, height_graph]);
 xy1.push([x_base+width_graph, height_graph]);
 var line = d3.svg.line()
       .x(function(d) {return d[0];})
       .y(function(d) {return d[1];});
 svg.append("path").attr("d", line(xy1) ).attr("fill","#F0E68C").attr("stroke","#D0C66C");
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");
 document.getElementById("optimal_Kernel").innerHTML = "Optimal bandwidth = <font color=\"red\">" + opt.toFixed(2) + "</font>";
 
 time_old[4] = new Date().getTime();
 
 // Kernel2(D)
 var wrap = d3.select('#graph_Kernel2');
 wrap.select("svg").remove(); // initialization
 var svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);

 var xy2 = new Array();
 for (var i = 0;i<opty2.length;i++) {
  xy2[i] = [x_base + Math.round(i*width_graph/(opty2.length-1)), height_graph - Math.round(height_graph*opty2[i]/(1.2*maxy))];
 }
 xy2.unshift([x_base, height_graph]);
 xy2.push([x_base+width_graph, height_graph]);
 var line = d3.svg.line()
       .x(function(d) {return d[0];})
       .y(function(d) {return d[1];});
 svg.append("path").attr("d", line(xy2) ).attr("fill","#FFDEAD").attr("stroke","#DFBE8D");
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");
 document.getElementById("optimal_Kernel2").innerHTML = "Optimal bandwidth = <font color=\"red\">" + opt.toFixed(2) + "</font>";
}

function DrawGraph_HMM(spike_time){
 var wrap = d3.select('#graph_HMM');
 wrap.select("svg").remove(); // initialization
 var svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);
 
 var x,y,maxy;
 var opty;
 var opt = ((offset-onset)/(spike_time.length-1)) * 5; // step width = ISI * 5
 opty = get_hmm_ratefunc(spike_time, opt);
 for(var i=0; i<opty.length; i++){
  if(i==0 || maxy<opty[i][1]) maxy=opty[i][1];
 }
 var x,y,xx,yy;
 for (var i = 0; i < opty.length; i++) {
  var x_pos=x_base+i*width_graph/opty.length;
  var height=height_hist*opty[i][1]/maxy;
     if (onset + i * opt < offset){
      svg.append("rect").attr("x", x_pos).attr("y", height_graph-height).attr("width", width_graph/opty.length+1).attr("height", height).attr("fill","#DA75F3");
     }else{
   x_pos=offset;
      svg.append("rect").attr("x", x_pos).attr("y", height_graph-height).attr("width", width_graph/opty.length+1).attr("height", height).attr("fill","#DA75F3");
     }
 }
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");
}

/*
DrawGraph_Bayes
estimates the firing rate with Kalman filtering and draw the rate. 

arguments:
spike_time: spike train

output
draw the estimated firing rates 

internal parameters
wrap, svg, maxy, xy: parameters for drawing a figure
kalman_data: firing rate estimated by Kalman filtering*/


function DrawGraph_Bayes(spike_time){
 var wrap = d3.select('#graph_Bayes');
 wrap.select("svg").remove(); // initialization
 var svg = wrap.append("svg").attr("width",x_base+width_graph).attr("height",height_graph);

 var maxy;
 var xy = new Array();
 
 var kalman_data = SecondStage(spike_time);
 // ThirdStage(spike_time,beta);
 for(var i=0; i<kalman_data[0].length; i++){
  if(i==0 || maxy<kalman_data[0][i]) maxy=kalman_data[0][i];
 }
 for (var i = 0;i<spike_time.length-1;i++) {
  xy[i] = [x_base + width_graph*(spike_time[i]/2+spike_time[i+1]/2-spike_time[0])/(spike_time[spike_time.length-1]-spike_time[0]), height_graph - height_graph*kalman_data[0][i]/(1.2*maxy)];
 }
 xy.unshift([x_base, height_graph - height_graph*kalman_data[0][0]/(1.2*maxy)]);
 xy.unshift([x_base, height_graph]);
 xy.push([x_base+width_graph, height_graph - height_graph*kalman_data[0][spike_time.length-2]/(1.2*maxy)]);
 xy.push([x_base+width_graph, height_graph]);
 var line = d3.svg.line()
       .x(function(d) {return d[0];})
       .y(function(d) {return d[1];});
 svg.append("path").attr("d", line(xy) ).attr("fill","#FFC0CB").attr("stroke","#DFA0AB");
 svg.append("rect").attr("x", x_base).attr("y", 0).attr("width", width_graph).attr("height", height_graph).attr("stroke","black").attr("stroke-width",1).attr("fill","none");
}
/* dft method */
function org_dft(x){
 var n = x.length;
 var y = new Array();// y[0] = y_re[]; y[1] = y_im[];
 y[0] = new Array();
 y[1] = new Array();
 for(var i=0;i<n;i++){
  y[0][i] = 0;
  y[1][i] = 0;
  for(var j=0;j<n;j++){
   y[0][i] += x[j]*Math.cos(2*Math.PI/n*i*j);
   y[1][i] += x[j]*(-Math.sin(2*Math.PI/n*i*j));
  }
 }
 return y;
}

function org_idft(y){
 var n = y[0].length;
 // input : y[0] = y_re[]; y[1] = y_im[];
 var w_re = Math.cos(-2*Math.PI/n);
 var w_im = -Math.sin(-2*Math.PI/n);
 var x = new Array();
 for(var i=0;i<n;i++){
  x[i] = 0;
  for(var j=0;j<n;j++){
   x[i]+= (y[0][j]*Math.cos(2*Math.PI/n*i*j) - y[1][j]*Math.sin(2*Math.PI/n*i*j))/n;
   // calculate real part only 
  }
 }
 return x;
}

function org_fft(n,re,im){
 var theta = 2*Math.PI/n;
 var firstr = new Array();
 var firsti = new Array();
 var secondr = new Array();
 var secondi = new Array();

 for(var i=0;i<n/2;i++){
  firstr[i] = re[2*i];
  firsti[i] = im[2*i];
 }
 for(var i=0;i<n/2;i++){
  secondr[i] = re[2*i+1];
  secondi[i] = im[2*i+1];
 }
 if(n/2>1){
  org_fft(n/2,firstr,firsti);
  org_fft(n/2,secondr,secondi);
 }
 for(var i=0;i<n/2;i++){
  var wr = Math.cos(theta * i);
  var wi = Math.sin(theta * i);
  re[i] = firstr[i] + wr*secondr[i] - wi*secondi[i];
  im[i] = firsti[i] + wr*secondi[i] + wi*secondr[i];

  wr = Math.cos(theta * (i+n/2));
  wi = Math.sin(theta * (i+n/2));
  re[i+n/2] = firstr[i] + wr*secondr[i] - wi*secondi[i];
  im[i+n/2] = firsti[i] + wr*secondi[i] + wi*secondr[i];
 }
}

function org_ifft(n,re,im){
 for(var i=0;i<n;i++){
  im[i] = -im[i];
 }
 org_fft(n,re,im);
 for(var i=0;i<n;i++){
  re[i] = re[i]/n
  im[i] = -im[i]/n;
 }
}

/*
Function kernel_rate returns optimized kernel densities estimate using a Gauss kernel function.

Input arguments
spike_time: sample data list.

Output arguments
y1: Estimated density using a Gauss kernel function.
y2: Estimated density using a Gauss kernel function with reflection boundary.
maxy: Maximum value of y2.
optw: Optimal kernel bandwidth.

Optimization principle: 
The optimal bandwidth is obtained as a minimizer of the fromula,
sum_{i, j} \int k(x - x_i) k(x - x_j) dx - 2 sum_{i~=j} k(x_i - x_j),
where k(x) is the kernel function, according to

Hideaki Shimazaki and Shigeru Shinomoto
Kernel Bandwidth Optimization in Spike Rate Estimation
Journal of Computational Neuroscience 2010
http://dx.doi.org/10.1007/s10827-009-0180-4

The above optimization is based on a principle of minimizing
expected L2 loss function between the kernel estimate and an
unknown underlying density function. An assumption is merely
that samples are drawn from the density independently each other.

For more information, please visit
http://2000.jukuin.keio.ac.jp/Shimazaki
*/

function kernel_rate(spike_time, y1, y2){

    var T = spike_time[spike_time.length-1] - spike_time[0];
    var dt_samp = spike_time[1]-spike_time[0];
    for (var i=0;i<spike_time.length-1;i++){
     if(dt_samp>spike_time[i+1]-spike_time[i]) dt_samp = spike_time[i+1]-spike_time[i];
    }
    var t_num=1000;
    if(Math.ceil(T/dt_samp)<t_num){
     t_num = Math.ceil(T/dt_samp);
    }
    var t = new Array();
    t[0] = spike_time[0]
    for(var i=0; i<t_num-1;i++){
     t[i+1]=t[i]+T/(t_num);
    }
 var dt = t[1]-t[0];
 for(var i=0;i<t.length-1;i++){
  if(dt>t[i+1]-t[i]){
   dt=t[i+1]-t[i];
  }
 }
    var y_hist = new Array();
    for(var i=0;i<t.length;i++){
     y_hist[i]=0;
    }
    for(var i=0;i<spike_time.length;i++){
     for(var j=0;j<t.length-1;j++){
         if(spike_time[i]>=t[j]-dt/2 && spike_time[i]<t[j+1]-dt/2) y_hist[j]++;
     }
     if(spike_time[i]>=t[t.length-1]-dt/2) y_hist[t.length-1]++;
    }
 var L = y_hist.length;
 var N = 0;
 for(var i=0;i<L;i++){
  N+=y_hist[i];
 }
 for (var i=0;i<t.length;i++){
  y_hist[i] = y_hist[i]/N/dt;   // density
 }

 var Wmin = 2*dt;
 var Wmax = 1*(spike_time[spike_time.length-1] - spike_time[0]);

 var tol = Math.pow(10,-5); 
 var phi = (Math.sqrt(5) + 1)/2;        //golden ratio

 // a = Wmin; b = Wmax;
 var a=ilogexp(Wmin);
 var b=ilogexp(Wmax);

 var c1 = (phi-1)*a + (2-phi)*b;
 var c2 = (2-phi)*a + (phi-1)*b;
 var f1 = kernel_cost_function(y_hist,N,logexp(c1),dt);
 var f2 = kernel_cost_function(y_hist,N,logexp(c2),dt);

 var k = 1;
 var W = new Array();
 var C = new Array();
 var optw;
 while (Math.abs(b-a) > tol*(Math.abs(c1)+Math.abs(c2)) && k <= 20){
  if (f1 < f2) {   
         b = c2;
         c2 = c1;

         c1 = (phi - 1)*a + (2 - phi)*b;
         
         f2 = f1;
         f1 = kernel_cost_function(y_hist,N,logexp(c1),dt);
         
         // 170926 fix
         // W[k] = Math.log(1+Math.exp(c1));
         W[k] = logexp(c1)
         
         C[k] = f1;
         //var optw = Math.log(1+Math.exp(c1));
         optw = logexp(c1);
         //y = yh1./sum(yh1.*dt);  //make the final output a density
  }else{
         a = c1;
         c1 = c2;
         
         c2 = (2 - phi)*a + (phi - 1)*b;
         
         f1 = f2;
         f2 = kernel_cost_function(y_hist,N,logexp(c2),dt);
                  //W[k] = Math.log(1+Math.exp(c2));
         W[k] = logexp(c2);
         C[k] = f2;
         
         //var optw = Math.log(1+Math.exp(c2));
         optw = logexp(c2);
         //y = yh2./sum(yh2.*dt);
  }
     k = k + 1;
 }
 var yh = new Array();
 yh = fftkernel(y_hist,optw/dt);
 var sum_yh = 0;
 for(var i=0;i<yh.length;i++){
  sum_yh += yh[i];
 }
 for(var i=0;i<yh.length;i++){
  y1[i] = yh[i]/sum_yh/dt;
 }
 /* reflection part */
 var yh2 = new Array();
 yh2 = fftkernel_ref(y_hist,optw/dt);
 sum_yh = 0;
 for(var i=0;i<yh2.length;i++){
  sum_yh += yh2[i];
 }
 var maxy = yh2[0]/sum_yh/dt;
 for(var i=0;i<yh2.length;i++){
  y2[i] = yh2[i]/sum_yh/dt;
  if(maxy<y2[i]) maxy = y2[i];
 }
 
 var res = new Array(2);
 res[0] = maxy;
 res[1] = optw;
 return res;
}

function kernel_cost_function(y_hist, N, w, dt){
 var yh = fftkernel(y_hist,w/dt);  // density

 var sumyh = 0;
 for(var i=0;i<yh.length;i++){
  sumyh += Math.pow(yh[i],2);
 }
 var sumyh_hist = 0;
 for(var i=0;i<yh.length;i++){
  sumyh_hist += yh[i]*y_hist[i];
 }
 
 // formula for density
 var C = sumyh*dt - 2* sumyh_hist*dt + 2*1/Math.sqrt(2*Math.PI)/w/N; 
 C = C * N* N;
 return C;
}

/*
fftkernel(x, width)

Function fftkernel applies the Gauss kernel smoother to an input signal using FFT algorithm.

Input argument
x: Sample signal vector
width: Kernel bandwidth (the standard deviation) in unit of the sampling resplution of x.

Output argument
re.slice(0, L): Smoothed signal.

MAY 5 / 23, 2012 Author Hideaki Shimazaki
RIKEN Brain Science Institute
http://2000.jukuin.keio.ac.jp/shimazaki
*/

function fftkernel(x, width){
 var L=x.length;
 var Lmax = Math.floor(L+3*width);
 var n=1;
 while(n<Lmax){
  n=n*2;
 }
 /*
 var x_buf=new Float64Array(n);
 for (var k=0;k<n;k++){
  x_buf[k]=0;
 }
 for (var k=0;k<x.length;k++){
  x_buf[k]=x[k];
 }
 var y_new = new Array();
 y_new = org_dft(x_buf);
 */
 var re=new Float64Array(n);
 var im=new Float64Array(n);
 for (var k=0;k<n;k++){
  re[k]=0;
  im[k]=0;
 }
 for (var k=0;k<x.length;k++){
  re[k]=x[k];
 }
 org_fft(n,re,im);
 
 var f_old = new Array();
 for (var k=0;k<n;k++){
  f_old[k]=k/n;
 }
 var f = new Array();
 var k=0;
 for (;k<Math.ceil(n/2)+1;k++){
  f[n-1-k]=f_old[k+1];
  f[k]=-f_old[k];
 }
 var K = new Array();
 for(var j=0;j<n;j++){
  K[j]=Math.exp(-0.5*Math.pow(width*2*Math.PI*f[j],2));
 }
 /*
 for(var j=0;j<n;j++){
  y_new[0][j] = y_new[0][j]*K[j];
  y_new[1][j] = y_new[1][j]*K[j];
 }
 var x_new = new Array();
 x_new = org_idft(y_new);
 return x_new.slice(0,L);
 */

 for(var j=0;j<n;j++){
  re[j] = re[j]*K[j];
  im[j] = im[j]*K[j];
 }
 org_ifft(n,re,im);
 return re.slice(0,L);
}

/*
fftkernel_ref(x, width)

Function fftkernel applies the Gauss kernel smoother to an input signal using FFT algorithm with reflection boundary.

Input argument
x: Sample signal vector
width: Kernel bandwidth (the standard deviation) in unit of the sampling resplution of x.

Output argument
y: Smoothed signal.

MAY 5 / 23, 2012 Author Hideaki Shimazaki
RIKEN Brain Science Institute
http://2000.jukuin.keio.ac.jp/shimazaki
*/

function fftkernel_ref(x,width){
 var yh = fftkernel(x,width);
 var halflen = Math.ceil(x.length/2);
 var remlen = x.length - halflen;
 var x_revleft = new Array();
 for(var i=0;i<remlen;i++){
  x_revleft[i] = 0;
 }
 for(var i=0;i<halflen;i++){
  x_revleft[remlen+i] = x[i];
 }
 var addleft = fftkernel(x_revleft,width);
 var x_revright = new Array();
 for(var i=0;i<halflen;i++){
  x_revright[i] = x[halflen+i];
 }
 for(var i=0;i<remlen;i++){
  x_revright[halflen+i] = 0;
 }
 var addright = fftkernel(x_revright,width);
 var y = new Array();
 for(var i=0;i<Math.ceil(yh.length/2);i++){
  y[i] = yh[i] + addleft[halflen-1-i];
  y[x.length-i-1] = yh[x.length-i-1] + addright[halflen+i];
 }
 return y;
}


function logexp(x){
 var y = new Array();
 if(x<100){
  y = Math.log(1+Math.exp(x));
 }else{
  y = x;
 }
 return y;
}

function ilogexp(x){
 // ilogexp = @(x) log(exp(x)-1);
 var y = new Array();
 if(x<100){
  y = Math.log(Math.exp(x)-1);
 }else{
  y = x;
 }
 return y;
}


function kern12(spike_time, width, y1, y2) {
 /*
 var x = new Array(res_graph)
 x[0] = onset;
 
 var maxy=0;
 var gauss;
 var addNumber = 0;
 
 for (var i=0; i<res_graph; i++) {
  x[i+1] = x[i] + (offset-onset)/(res_graph-1); 
 }
 var maxy=0;
 var gauss;
 for (var i=0; i<res_graph; i++) {
  addNumber = 0;
  y1[i] = 0;
  for (var j in spike_time) {
   if((x[i]-5*width <= spike_time[j]) && (spike_time[j] <= x[i]+5*width)){
    gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-spike_time[j])*(x[i]-spike_time[j])/2/width/width);
    y1[i] = y1[i] + gauss / spike_time.length;
   }
   if (x[i]-5*width < onset) {
    if (-(x[i]-5*width)+2*onset > spike_time[j]){
     gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-(onset-(spike_time[j]-onset)))*(x[i]-(onset-(spike_time[j]-onset)))/2/width/width);
     addNumber = addNumber + gauss / spike_time.length;
    }
   }else if(x[i]+5*width > offset){
    if(-(x[i]+5*width)+2*offset > spike_time[i]){
     gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-(offset+(offset-spike_time[j])))*(x[i]-(offset+(offset-spike_time[j])))/2/width/width);
     addNumber = addNumber + gauss / spike_time.length;
    }
   }
  }
  y2[i] = y1[i] + addNumber;
  if(maxy<y2[i]) maxy=y2[i];
 }
*/
 //////////////////
  var L=spike_time.length;
  var Lmax = L+3*width;
  var n=1;
  while(n<Lmax){
   n=n*2;
  }
  var imag=new Float64Array(n);
  var real=new Float64Array(n);
  for (var k=0;k<n;k++){
   imag[k]=0;
   real[k]=0;
  }
  for (var k=0;k<spike_time.length;k++){
   real[k]=spike_time[k];
  }
  fftnoasm=new FftModule(n,false);
  fftnoasm.fft(real,imag,0);
  var f_old = new Array();
  for (var k=0;k<real.length;k++){
   f_old[k]=k/real.length;
  }
  var f = new Array();
  var k=0;
  for (;k<Math.ceil(real.length/2)+1;k++){
   f[real.length-1-k]=f_old[k+1];
   f[k]=-f_old[k];
  }
  var K = new Array();
  for(var j=0;j<real.length;j++){
   K[j]=Math.exp(-0.5*Math.pow(width*2*Math.PI*f[j],2));
  }
  for(var j=0;j<real.length;j++){
   y1[j] = real[j]*K[j];
   imag[j] = imag[j]*K[j];
  }
  fftnoasm.fft(y1,imag,1);
  document.data.spikes.value = y1;

  var maxy = 0;
  for(var i = 0;i<y1.length;i++){
   if(maxy<y1[i]) maxy=y1[i];
  }
 

 //////////////////
 return maxy;
}

// estimate the firing rate, given the parameters of binsize and binrate.
function EstimateRate(spike_time, opt_binsize, opt_rate) {
  var opt_binnum = Math.ceil((spike_time[spike_num - 1] - onset) / opt_binsize);
  var rate_max;

  for (var i = 0; i < opt_binnum; i++) {
    opt_rate[i] = 0;
  }
  for (i = 0; i < spike_num; i++) {
    opt_rate[Math.floor((spike_time[i] - onset) / opt_binsize)] += 1.0 / opt_binsize;
  }
  for (i = 0; i < opt_binnum; i++) {
    if (i == 0 || opt_rate[i] > rate_max) rate_max = opt_rate[i];
  }
  return rate_max;
}
// output
function GenerateOutputFileMessage(message) {
 return "<div id='Output'></div> <script type='text/javascript'>var myBlob = new Blob([\"" + message + "\"], {type: 'text/html'}); var url = URL.createObjectURL(myBlob); document.getElementById('Output').innerHTML = '<a href=' + url + ' download=datasheet.csv>download as csv</a>';</script>";
}

function OutputResults_SS() {
 var result;
 var spike_time = new Array();
 PostData(spike_time);
 var opt_binsize = new Array();
 var opt_rate = new Array();
 opt_binsize = SSOS(spike_time);
 EstimateRate(spike_time, opt_binsize[0], opt_rate);

 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += onset.toFixed(2) + ",0\\n";
 for (var i = 0; i < opt_rate.length; i++) {
  filemessage += (onset + i * opt_binsize[0]).toFixed(2) + "," + opt_rate[i].toFixed(2) + "\\n";
  filemessage += (onset + (i + 1) * opt_binsize[0]).toFixed(2) + "," + opt_rate[i].toFixed(2) + "\\n";
 }
 filemessage += (onset + opt_rate.length * opt_binsize[0]).toFixed(2) + ",0\\n";

 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Poissonian optimization</h2>");
 WIN_RESULTS.document.writeln("Optimal binsize: <b>"+opt_binsize[0].toFixed(2)+"</b><br><br>");

 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage));

 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td width=150> Y-AXIS (density)</td>");
 WIN_RESULTS.document.writeln("<tr align=right><td>" + onset.toFixed(2) + "</td><td>0.00</td></tr>");
 for (var i=0;i<opt_rate.length;i++) {
  WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + i * opt_binsize[0]).toFixed(2) + "</td><td>" + opt_rate[i].toFixed(2) + "</td></tr>");
  WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + (i + 1) * opt_binsize[0]).toFixed(2) + "</td><td>" + opt_rate[i].toFixed(2) + "</td></tr>");
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + opt_rate.length * opt_binsize[0]).toFixed(2) + "</td><td>0.00</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.close();
}

function OutputResults_OS() {
 var result;
 var spike_time = new Array();
 PostData(spike_time);
 var opt_binsize = new Array();
 var opt_rate = new Array();
 opt_binsize = SSOS(spike_time);
 EstimateRate(spike_time, opt_binsize[1], opt_rate);

 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += onset.toFixed(2) + ",0\\n";
 for (var i = 0; i < opt_rate.length; i++) {
  filemessage += (onset + i * opt_binsize[1]).toFixed(2) + "," + opt_rate[i].toFixed(2) + "\\n";
  filemessage += (onset + (i + 1) * opt_binsize[1]).toFixed(2) + "," + opt_rate[i].toFixed(2) + "\\n";
 }
 filemessage += (onset + opt_rate.length * opt_binsize[1]).toFixed(2) + ",0\\n";

 
 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Non-Poissonian optimization</h2>");
 WIN_RESULTS.document.writeln("Optimal binsize: <b>"+opt_binsize[1].toFixed(2)+"</b><br>Lv: <b>" + lv.toFixed(2) + "</b> (" + np + " firing)<br><br>");

 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage)); 

 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td width=150> Y-AXIS (density) </td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>" + onset.toFixed(2) + "</td><td>0.00</td></tr>");
 for (var i=0;i<opt_rate.length;i++) {
  WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + i * opt_binsize[1]).toFixed(2) + "</td><td>" + opt_rate[i].toFixed(2) + "</td></tr>");
  WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + (i + 1) * opt_binsize[1]).toFixed(2) + "</td><td>" + opt_rate[i].toFixed(2) + "</td></tr>");
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>" + (onset + opt_rate.length * opt_binsize[1]).toFixed(2) + "</td><td>0.00</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.close();
}

function xaxisForKernel(spike_time) {
 var x = new Array(res_graph);
 var data_max = spike_time[spike_time.length - 1];
 var data_min = spike_time[0];
 x[0] = data_min;
 for (var i = 0; i < res_graph - 1; i++) {
  x[i + 1] = x[i] + (data_max - data_min) / (res_graph - 1);
 }
 return x;
}

function kern(spike_time, width, y) {
 var x = new Array(res_graph)
 x[0] = onset;
 for (var i=0; i<res_graph; i++) {
  x[i+1] = x[i] + (offset-onset)/(res_graph-1); 
 }
 var maxy=0;
 var gauss;
 for (var i=0; i<res_graph; i++) {
  y[i] = 0;
  for (var j in spike_time) {
   if((x[i]-5*width <= spike_time[j]) && (spike_time[j] <= x[i]+5*width)){
    gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-spike_time[j])*(x[i]-spike_time[j])/2/width/width);
    y[i] = y[i] + gauss;
   }
  }
  if(maxy<y[i]) maxy=y[i];
 }
 return maxy;
}

function kern2(spike_time, width, y) {
 var x = new Array(res_graph)
 x[0] = onset;
 for (var i=0; i<res_graph; i++) {
  x[i+1] = x[i] + (offset-onset)/(res_graph-1); 
 }
 var maxy=0;
 var gauss;
 var addNumber = 0;

 for (var i=0; i<res_graph; i++) {
  addNumber = 0;
  y[i] = 0;
  for (var j in spike_time) {
   if((x[i]-5*width <= spike_time[j]) && (spike_time[j] <= x[i]+5*width)){
    gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-spike_time[j])*(x[i]-spike_time[j])/2/width/width);
    y[i] = y[i] + gauss;
   }
  
    
   if (x[i] - 5*width<onset) {
    if (-(x[i]-5*width)+2*onset > spike_time[j]){
     gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-(onset-(spike_time[j]-onset)))*(x[i]-(onset-(spike_time[j]-onset)))/2/width/width);
     addNumber = addNumber + gauss;
    }
   }else if(x[i]+5*width>offset){
    if(-(x[i]+5*width)+2*offset > spike_time[i]){
     gauss = 1/Math.sqrt(2*Math.PI)/width*Math.exp(-(x[i]-(offset+(offset-spike_time[j])))*(x[i]-(offset+(offset-spike_time[j])))/2/width/width);
     addNumber = addNumber + gauss;
    }
   }
  }
  y[i] += addNumber;
  if(maxy<y[i]) maxy=y[i];
 }
 return maxy;
}


function OutputResults_Kernel() {
 var spike_time = new Array();
 PostData(spike_time);
 var opt = Kernel(spike_time);
 var opty = new Array();
 kern(spike_time, opt, opty);
 var xaxis = xaxisForKernel(spike_time);

 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += xaxis[0].toFixed(3) + ",0\\n";
 for (var i = 0; i < xaxis.length; i++) {
  filemessage += xaxis[i].toFixed(3) + "," + opty[i].toFixed(3) + "\\n";
 }
 filemessage += xaxis[xaxis.length - 1].toFixed(3) + ",0\\n";

 
 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Kernel Density Estimation</h2>");
 WIN_RESULTS.document.writeln("Optimal Bandwidth: <b>"+opt.toFixed(3)+"</b><br><br>");

 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage));
 
 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td width=150> Y-AXIS (density) </td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[0].toFixed(3)+"</td><td>0.00</td></tr>");
 for (var i=0;i<xaxis.length;i++) {
  WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[i].toFixed(3)+"</td><td>" + opty[i].toFixed(3) + "</td></tr>");
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[xaxis.length-1].toFixed(3)+"</td><td>0.00</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.close();
}

function OutputResults_Kernel2() {
 var spike_time = new Array();
 PostData(spike_time);
 var opt = Kernel(spike_time);
 var opty = new Array();
 kern2(spike_time, opt, opty);
 var xaxis = xaxisForKernel(spike_time);

 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += xaxis[0].toFixed(3) + ",0\\n";
 for (var i = 0; i < xaxis.length; i++) {
  filemessage += xaxis[i].toFixed(3) + "," + opty[i].toFixed(3) + "\\n";
 }
 filemessage += spike_time[spike_time.length - 1] + ",0\\n";
 
 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Kernel Density Estimation</h2>");
 WIN_RESULTS.document.writeln("Optimal Bandwidth: <b>"+opt.toFixed(3)+"</b><br><br>");
 
 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage));
 
 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td> Y-AXIS (density) </td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[0].toFixed(3)+"</td><td>0.00</td></tr>");
 for (var i=0;i<xaxis.length;i++) {
  WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[i].toFixed(3)+"</td><td>" + opty[i].toFixed(3) + "</td></tr>");
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>"+xaxis[xaxis.length -1].toFixed(3)+"</td><td>0.00</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.close();
}

function OutputResults_HMM() {
 var spike_time = new Array();
 PostData(spike_time);
 var opty;
 var opt = (offset-onset)/(spike_time.length-1);
 var time = onset;
 opty = get_hmm_ratefunc(spike_time, opt); // step size = 5* (mean inter-spike interval)
 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += time.toFixed(3) + ",0\\n";
 filemessage += time.toFixed(3) + "," + opty[0][1].toFixed(3) + "\\n";
 time += opt;
 for (var i = 1; i < opty.length; i++) {
  if (opty[i][1] != opty[i - 1][1]) {
   filemessage += time.toFixed(3) + "," + opty[i - 1][1].toFixed(3) + "\\n";
   filemessage += time.toFixed(3) + "," + opty[i][1].toFixed(3) + "\\n";
  }
  time += opt;
 }
 filemessage += time.toFixed(3) + "," + opty[opty.length - 1][1].toFixed(3) + "\\n";
 filemessage += time.toFixed(3) + ",0\\n";
 time = onset;
     
 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Two state hidden Markov model</h2>");

 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage));

 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td width=150> Y-AXIS (density) </td></tr>");
 //WIN_RESULTS.document.writeln("<tr align=right><td>0.000</td><td>0.000</td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>" +time.toFixed(2)+"</td><td>0.000</td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>" + time.toFixed(2) + "</td><td>" + opty[0][1].toFixed(3) + "</td></tr>");
 time+=opt;
 for (var i=1;i<opty.length;i++) {
  if(opty[i][1]!=opty[i-1][1]){
   WIN_RESULTS.document.writeln("<tr align=right><td>"+time.toFixed(2)+"</td><td>" + opty[i-1][1].toFixed(3)+"</td></tr>");
   WIN_RESULTS.document.writeln("<tr align=right><td>"+time.toFixed(2)+"</td><td>" + opty[i][1].toFixed(3) + "</td></tr>");
  }
  time+=opt;
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>"+ time.toFixed(2) +"</td><td>" + opty[opty.length-1][1].toFixed(3) + "</td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>"+ time.toFixed(2) +"</td><td>0.000</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.writeln("</blockquote>");
 WIN_RESULTS.document.close();
}

function OutputResults_Bayes(){
 var spike_time = new Array();
 PostData(spike_time);
 var opty;
 var kalman_data = SecondStage(spike_time);
 // ThirdStage(spike_time,beta);

 //save as csv
 var filemessage = "X-AXIS,Y-AXIS\\n";
 filemessage += ((spike_time[0] + spike_time[1]) / 2).toFixed(3) + ",0\\n";
 for (var i = 0; i < spike_time.length - 1; i++) {
  filemessage += ((spike_time[i] + spike_time[i + 1]) / 2).toFixed(3) + "," + kalman_data[0][i].toFixed(3) + "\\n";
 }
 filemessage += ((spike_time[spike_time.length - 2] + spike_time[spike_time.length - 1]) / 2).toFixed(3) + ",0\\n";

 WIN_RESULTS = window.open();
 WIN_RESULTS.document.open();
 WIN_RESULTS.document.writeln("<title>Data Sheet of the Optimized Histogram</title>");
 WIN_RESULTS.document.writeln("<h2>Histgram: Bayesian model Estimation</h2>");
 WIN_RESULTS.document.writeln("<br><br>");

 WIN_RESULTS.document.writeln(GenerateOutputFileMessage(filemessage));
 
 WIN_RESULTS.document.writeln("<table border=1><tr align=center><td width=150> X-AXIS (time)  </td><td width=150> Y-AXIS (density) </td></tr>");
 WIN_RESULTS.document.writeln("<tr align=right><td>"+((spike_time[0] + spike_time[1]) / 2).toFixed(3)+"</td><td>0.00</td></tr>");
 for (var i=0;i<spike_time.length - 1;i++) {
  WIN_RESULTS.document.writeln("<tr align=right><td>"+((spike_time[i] + spike_time[i + 1]) / 2).toFixed(3)+"</td><td>" + kalman_data[0][i].toFixed(3) + "</td></tr>");
 }
 WIN_RESULTS.document.writeln("<tr align=right><td>"+((spike_time[spike_time.length - 2] + spike_time[spike_time.length - 1]) / 2).toFixed(3) + "</td><td>0.00</td></tr>");
 WIN_RESULTS.document.writeln("</table><br>");
 WIN_RESULTS.document.close();
}
