clc;
clear all;

branches = [%From	To	 R+jX (ohms)
1	2	0
2	3	0.0232+0.0475i
3	4	0.0441+0.0447i
3	5	0.0629+0.0638i
3	9	0.0174+0.0356i
5	6	0.0503+0.051i
5	7	0.0818+0.0829i
7	8	0.0629+0.0638i
9	10	0.0116+0.0238i
10	14	0.0566+0.0574i
10	11	0.0566+0.0574i
10	15	0.0174+0.0356i
11	126	0
126	16	0.1069+0.10846i
15	36	0.0378+0.0383i
15	20	0.0478+0.098i
16	13	0.0629+0.0638i
16	12	0.0629+0.0638i
17	18	0.0944+0.0957i
17	19	0.0881+0.0893i
20	21	0.0629+0.0638i
20	23	0.0174+0.0356i
21	22	0.0818+0.0829i
23	24	0.1322+0.134i
23	25	0.0145+0.0297i
25	26	0.1385+0.1404i
25	27	0.0159+0.0327i
27	127	0
127	28	0.0203+0.0455i
27	30	0.0116+0.0238i
28	29	0.0159+0.0358i
28	33	0.0566+0.0574i
29	35	0.1259+0.1276i
30	31	0.0174+0.0356i
31	32	0.0203+0.0416i
32	119	0.0116+0.0238i
33	34	0.0755+0.0766i
36	17	0.0252+0.0255i
37	38	0.0377+0.0845i
37	42	0.0145+0.0297i
38	39	0.0755+0.0766i
38	40	0.0629+0.0638i
40	41	0.0818+0.0829i
42	43	0.0818+0.0829i
42	44	0.0145+0.0297i
44	45	0.1259+0.1276i
44	46	0.0116+0.0238i
46	47	0.0503+0.051i
46	49	0.0145+0.0297i
47	48	0.0755+0.0766i
49	50	0.0087+0.0178i
49	51	0.0145+0.0297i
51	52	0.0145+0.0297i
52	53	0.0145+0.0297i
54	55	0.0116+0.0238i
55	56	0.0072+0.0148i
56	57	0.0159+0.0327i
56	59	0.0203+0.0416i
57	58	0.0159+0.0327i
59	60	0.0629+0.0638i
59	62	0.0435+0.0891i
60	61	0.0629+0.0638i
62	63	0.0319+0.0653i
62	64	0.048+0.0229i
64	65	0.0336+0.016i
65	66	0.0672+0.032i
66	67	0.0816+0.0389i
67	68	0.0624+0.0297i
69	70	0.0503+0.051i
69	74	0.0159+0.0327i
69	99	0.0145+0.0297i
70	71	0.0692+0.0702i
71	72	0.0818+0.0829i
72	73	0.0692+0.0702i
74	75	0.0692+0.0702i
74	78	0.0116+0.0238i
75	76	0.0881+0.0893i
76	77	0.1007+0.1021i
78	79	0.0232+0.0475i
78	88	0.0406+0.0831i
79	80	0.0058+0.0119i
80	81	0.013+0.0267i
80	82	0.0275+0.0564i
82	83	0.0101+0.0208i
83	84	0.0145+0.0297i
83	86	0.1699+0.1723i
84	85	0.0145+0.0297i
86	87	0.1196+0.1212i
88	89	0.0261+0.0534i
89	90	0.0441+0.0447i
89	91	0.0159+0.0327i
91	92	0.0566+0.0574i
91	93	0.013+0.0267i
93	94	0.0755+0.0766i
93	95	0.013+0.0267i
95	96	0.0692+0.0702i
95	97	0.0174+0.0356i
97	98	0.0503+0.051i
99	100	0.0159+0.0327i
100	101	0.0319+0.0653i
101	102	0.0174+0.0356i
102	123	0.0464+0.095i
103	104	0.0566+0.0574i
103	107	0.0159+0.0327i
104	105	0.0818+0.0829i
105	106	0.1762+0.1786i
107	108	0.0566+0.0574i
107	110	0.0188+0.0386i
108	109	0.1448+0.1467i
110	111	0.1133+0.1148i
110	121	0.05797+0.11876i
111	112	0.0755+0.0766i
112	113	0.1448+0.1467i
112	114	0.0315+0.0319i
114	115	0.1322+0.134i
115	116	0.0818+0.0829i
118	37	0.0217+0.0445i
117	54	0.0232+0.0475i
53	120	0.029+0.0594i
125	128	0
128	69	0.0203+0.0416i
124	103	0.0145+0.0297i
15	117	1e-06+1e-06i
20	118	1e-06+1e-06i
62	125	1e-06+1e-06i
99	124	1e-06+1e-06i
63	122	1e-05+0.0001i
    ];

spotloads = [%Bus Number	P+jQ(KVA)
1	0
2	0
3	13.333+6.6667i
4	6.6667+3.3333i
5	0
6	13.333+6.6667i
7	6.6667+3.3333i
8	13.333+6.6667i
9	6.6667+3.3333i
10	0
11	13.333+6.6667i
12	6.6667+3.3333i
13	13.333+6.6667i
14	6.6667+3.3333i
15	0
16	0
17	0
18	13.333+6.6667i
19	6.6667+3.3333i
20	0
21	13.333+6.6667i
22	13.333+6.6667i
23	0
24	13.333+6.6667i
25	0
26	13.333+6.6667i
27	0
28	0
29	0
30	13.333+6.6667i
31	13.333+6.6667i
32	13.333+6.6667i
33	6.6667+3.3333i
34	6.6667+3.3333i
35	13.333+6.6667i
36	13.333+6.6667i
37	13.333+6.6667i
38	0
39	13.333+6.6667i
40	6.6667+3.3333i
41	6.6667+3.3333i
42	0
43	6.6667+3.3333i
44	6.6667+3.3333i
45	13.333+6.6667i
46	0
47	6.6667+3.3333i
48	6.6667+3.3333i
49	35+25i
50	70+50i
51	46.667+31.667i
52	13.333+6.6667i
53	6.6667+3.3333i
54	13.333+6.6667i
55	13.333+6.6667i
56	0
57	6.6667+3.3333i
58	6.6667+3.3333i
59	0
60	6.6667+3.3333i
61	6.6667+3.3333i
62	6.6667+3.3333i
63	0
64	13.333+6.6667i
65	13.333+6.6667i
66	25+11.667i
67	46.667+33.334i
68	25+11.667i
69	0
70	6.6667+3.3333i
71	13.333+6.6667i
72	6.6667+3.3333i
73	13.333+6.6667i
74	0
75	13.333+6.6667i
76	13.333+6.6667i
77	13.333+6.6667i
78	81.666+60.001i
79	13.333+6.6667i
80	0
81	13.333+6.6667i
82	13.333+6.6667i
83	0
84	13.333+6.6667i
85	6.6667+3.3333i
86	6.6667+3.3333i
87	13.333+6.6667i
88	6.6667+3.3333i
89	13.333+6.6667i
90	13.333+6.6667i
91	0
92	13.333+6.6667i
93	0
94	13.333+6.6667i
95	0
96	13.333+6.6667i
97	6.6667+3.3333i
98	6.6667+3.3333i
99	0
100	13.333+6.6667i
101	13.333+6.6667i
102	13.333+6.6667i
103	0
104	6.6667+3.3333i
105	13.333+6.6667i
106	13.333+6.6667i
107	0
108	13.333+6.6667i
109	13.333+6.6667i
110	0
111	13.333+6.6667i
112	0
113	6.6667+3.3333i
114	6.6667+3.3333i
115	13.333+6.6667i
116	6.6667+3.3333i
117	0
118	0
119	0
120	0
121	0
122	0
123	0
124	0
125	0
126	0
127	0
128	0];

spotloads=[spotloads(:,1) spotloads(:,2)*1e3];%P+jQ(VA)
Vb = 4.16/sqrt(3);%(kV)
Vb = Vb*1e3;%(V)
rn = 1;% root node
[V, I] = FBS_LF(branches ,rn ,Vb,spotloads);
voltage = abs(V(:,2)/Vb);
plot(abs(V(:,1)),abs(V(:,2))/Vb)
%set(get(gca,'YAxis'),'Exponent',3)
ylabel('Voltage (V)')
xlabel('bus number')

for i=1:128
S(i,1) = V(i,2)*conj(I(i,2));
end

RP = real(S);
QP = imag(S);