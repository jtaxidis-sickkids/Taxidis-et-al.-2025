clear asap

%% SOM-CRE MICE
asap(1).name = 'MM7-1';     asap(1).ASAP3flag = 1;
asap(2).name = 'MM7-2';     asap(2).ASAP3flag = 1;
asap(3).name = 'ASAP3-1';   asap(3).ASAP3flag = 1;
asap(4).name = 'ASAP3-3';   asap(4).ASAP3flag = 1;
asap(5).name = 'ASAP3-4';   asap(4).ASAP3flag = 1;


%% SESSION DATA
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%                       DAY         VIDEOS TO POOL          VIDEOS TO KEEP     MAX TRIAL            ROI                   FOV               MINFREQ
% ---------------------------------------------------------------------------------------------------------------------------------------------------
asap(1).sessions = {'08_27_2020',   [],                     [2],                [8],                {},                   [],               [];          
                    '09_04_2020',   [],                     [1 3],              [8 6],              {[1],[]},             [],               [];              
                    '09_21_2020',   [2 5],                  [2 6],              [32 8],             {},                   [],               [];
                    '09_22_2020',   [1 2; 3 4],             [1 3 6],            [16 16 8],          {},                   [nan 1],          [1/3, 3, 1/3];
                    '09_23_2020',   [],                     [1 3],              [8 8],              {},                   [1 nan],          [5, 5]};   

                    
asap(2).sessions = {'08_28_2020',   [],                     [1 3 4],            [8 8 8],            {},                   [],               [1/3, 1/3, 1];               
                    '09_04_2020',   [1 2],                  [1 3],              [20 6],             {},                   [],               [1, 1/3, 1/3];               
                    '09_22_2020',   [2 3],                  [2],                [16],               {},                   [],               [2]}; 


asap(3).sessions = {'04_21_2021',   [1 4; 5 7],             [1 5],              [32 23],            {},                   [],               [2, 1/3];               
                    '04_24_2021',   [1 2; 3 4; 5 6],        [1 3 5],            [20 20 20],         {[],[],[],[],[1]},    [],               [];              
                    '04_26_2021',   [3 4],                  [1 2 3],            [8 8 16],           {},                   [],               [];  
                    '05_10_2021',   [2 4; 5 6],             [1 2 5],            [8 28 24],          {[1],[],[]},          [1 2 3],          [1/3, 1/3, 3];  
                    '05_11_2021',   [5 7],                  [1 4 5],            [12 10 36],         {},                   [1 nan nan],      [];  
                    '05_13_2021',   [],                     [3 4 6],            [12 12 20],         {},                   [3 4 5],          [];  
                    '05_18_2021',   [],                     [1 2 5 6 8],        [8 12 8 12 8],      {[],[],[3],[],[]},    [1 2 3 4 5],      [4, 1/3, 1, 1/3, 1/3]}; 
 
                
asap(4).sessions = {'04_22_2021',   [1 3; 4 5; 6 8],        [1 4 6 9 10],       [24 16 24 8 8],     {[] [2],[],[],[]},    [],               [];               
                    '04_24_2021',   [5 8; 9 10; 13 14],     [1 3 5 9 11 13 15], [10 8 32 16 8 16 10], {[],[],[],[],[],[1],[]},  [],         [];               
                    '05_08_2021',   [1 4; 7 8; 9 11],       [1 5 7 9],          [48 12 24 36],      {},                   [1 2 3 4],        [];  % In sess7-8 3 cells all spiking!!!!
                    '05_10_2021',   [5 6; 7 8],             [2 3 5 7],          [12 12 24 24],      {},                   [],               [];  
                    '05_11_2021',   [2 3],                  [1 2 5 6],          [10 40 12 12],      {},                   [nan 1 2 nan],    [];  
                    '05_14_2021',   [2 3],                  [1 2 4],            [12 20 12],         {},                   [nan 2 3],        [1/3, 2, 1/3]; 
                    '05_15_2021',   [],                     [1 3 5 6 8],        [12 12 12 8 8],     {},                   [1 2 3 4 5],      [];  
                    '05_18_2021',   [],                     [1 4 6 9],          [8 8 8 8],          {},                   [1 2 3 4],        [1/3, 1, 1/3, 1.3]}; 
 
asap(5).sessions = {'03_30_2022',   [9 10],                 [1 2 4 5 7 8 9],    [12 12 12 10 8 10 28], {},                [1 2 3 4 6 7 8],  [];
                    '04_01_2022',   [2 3],                  [1 2 4],            [12 20 12],         {},                   [1 2 3],          []; 
                    '04_13_2022',   [8 9],                  [1 2 3 8],          [10 12 11 20],      {},                   [1 2 3 nan],      []; 
                    '04_14_2022',   [],                     [1],                [10],               {},                   [1],              []; 
                    '04_15_2022',   [],                     [3],                [8],                {},                   [2],              [] }; 
 
%% SAME CELLS
% --------------------------------------------------------------------------------------------------------------------------------------
%                      DAY   SESSION   ROI
% --------------------------------------------------------------------------------------------------------------------------------------
% MM7-1 -------------------------------------------------------------------
asap(1).samecells{1} = [4,      3,      1;
                        5,      1,      1];
% -------------------------------------------------------------------------

% ASAP3-1 -----------------------------------------------------------------
asap(3).samecells{1} = [4,      5,      1;
                        5,      1,      1;
                        6,      3,      1;
                        7,      8,      1]; 
                    
asap(3).samecells{2} = [6,      4,      1;
                        7,      1,      1]; 

asap(3).samecells{3} = [6,      6,      1;
                        7,      2,      1];  
% -------------------------------------------------------------------------

% ASAP3-3 -----------------------------------------------------------------
asap(4).samecells{1} = [3,      1,      1;
                        5,      2,      1];

asap(4).samecells{2} = [6,      2,      1;
                        7,      5,      1];

asap(4).samecells{3} = [6,      4,      1;
                        7,      1,      1;
                        8,      1,      1];

asap(4).samecells{4} = [7,      3,      1;
                        8,      4,      1];
% -------------------------------------------------------------------------

% ASAP3-4 -----------------------------------------------------------------
asap(5).samecells{1} = [1,      8,      1;
                        3,      2,      1];

asap(5).samecells{2} = [2,      2,      1;
                        3,      3,      1];

asap(5).samecells{3} = [2,      1,      1;
                        4,      1,      1];
% -------------------------------------------------------------------------






