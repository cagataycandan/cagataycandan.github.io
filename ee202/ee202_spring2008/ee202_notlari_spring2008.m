a=[
76	99	85
83	89	83
80	83	87
74	73	91
88	82	72
70	69	84
58	80	77
55	77	77
76	67	67
46	75	81
73	70	60
77	65	56
65	56	70
56	46	85
56	62	70
64	50	71
66	65	55
67	54	60
51	54	69
62	64	38
64	36	60
49	47	61
55	57	46
61	41	50
63	33	56
54	60	33
48	19	70
16	40	77
18	43	66
57	43	36
51	36	42
57	44	27
64	22	40
54	24	39
23	48	42
37	30	47
45	20	46
39	51	23
34	43	26
56	21	28
53	17	32
36	17	43
38	33	24
46	17	33
26	16	45
52	38	0
28	18	1
19	13	1
19	0	0
0	2	8
4	0	0];    

numstu=size(a,1);
a2=a-repmat(mean(a),numstu,1);
q=a2'*a2/numstu;d=sqrt(diag(q));a3=a2*inv(diag(d));
a3'*a3/numstu

figure(1),plot(a(:,1),a(:,3),'x');xlabel('MT1');ylabel('Final'); %print -djpeg90 -f1 mt1_final.jpeg
figure(2),plot(a(:,1),a(:,2),'x');xlabel('MT1');ylabel('MT2'); %print -djpeg90 -f1 mt1_mt2.jpeg
figure(3),plot(a(:,2),a(:,3),'x');xlabel('MT2');ylabel('Final'); %print -djpeg90 -f1 mt2_final.jpeg
