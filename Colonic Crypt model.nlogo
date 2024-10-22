extensions [csv]
GLOBALS
[
  hours                     ;TIME UNIT OF THE MODEL 1 HOUR = 1 TICK;
  cells-line                ;LINE OF ALL CELLS IN THAT LINE;
  max-cells-line            ;THE MAXIMUM CELLS' LINE OF ALL CELLS;
  n-cells-dead              ;NUMBER OF DEAD CELLS FOR CONDITION a) AND g)
  n-hypoxia-death              ;NUMBER OF DEAD CELLS FOR CONDITION d) AND e)
  met-death                 ;NUMBER OF DEAD CELLS FOR CONDITION f)
  immune-death              ;NUMBER OF DEAD CELLS FOR CONDITION e)
  energy-death              ;NUMBER OF DEAD CELLS FOR CONDITION h) AND i)
  n-tumoral-dead            ;NUMBER OF DEAD CELLS FOR CONDITION c)
  n-repetitions             ;NUMBER OF PROCEDURES REPEATED BY CELLS FOR EVERY LINE;
  APC-mutations             ;COUNT THE NUMBER OF MUTATION (ITEM CHANGE) OF APC GENE
  Kras-mutations            ;COUNT THE NUMBER OF MUTATION (ITEM CHANGE) OF KRAS GENE
  p53-mutations             ;COUNT THE NUMBER OF MUTATION (ITEM CHANGE) OF P53GENE
  number-other-mut          ;COUNT THE SUM OF MUTATION (ITEM CHANGE) OF GENES GENE
  crypt
  cells-dead
  count-cells
  p53-dead
]

breed [cells cell]
breed [immune-cells immune-cell]

patches-own
[
  hypoxia                     ;GRADIENT OF OXIGEN, DEPENDS ON THE NUMBER OF CELLS;
]

immune-cells-own
[
  task                      ;COUNT THE NUMBER OF CELLS KILLED BY IMMUNE SYSTEM CELLS;
  immune-life
]

cells-own
[
  tipology                  ;THERE ARE TWO TYPOLOGIES "TRANSIT-AMPLIFYING" AND "DIFFERENTIATED";
  life-time                 ;LIFESPAN OF THE CELLS;
  mitosis-time              ;MOMENT OF MITOTIC ACTIVITY;
  mitosis-mean-time         ;MEAN VALUE OF THE DISTRIBUTION OF MITOTICS TIME;
  B-cat                     ;GRADIENT IN FUNCTION OF WHICH DEPENDS THE DIFFERENTIATION OF CELLS FROM "TRANSIT-AMPLIFYING" TO "DIFFERENTIATED";
  B-block?                  ;IF APC IS MUTATED CHANGE THE STATE IN "TRUE" AND PREVENT THE CELLS DIFFERENTIATION;
  own-line                  ;EVERY CELL COUNT ITS OWN LINE;
  pre-adenoma?              ;BECOME "TRUE" IF APC IS MUTATED;
  adenoma?                  ;BECOME "TRUE" IF APC AND KRAS ARE MUTATED;
  tumoral?                  ;BECOME "TRUE" IF APC AND KRAS AND P53 ARE MUTATED;
  mitosis?                  ;BECOME "TRUE" IF LIFE-TIME IS A MULTIPLE OF THE MITOTIS-TIME;
  p53-control               ;IF P53 IS WILD-TYPE KILLS THE CELLS WITH SUM-GENES > THRESHOLD;
  APC                       ;LIST THAT REPRESENT APC GENE;
  Kras                      ;LIST THAT REPRESENT KRAS GENE;
  p53                       ;LIST THAT REPRESENT P53 GENE;
  genes                     ;LIST THAT REPRESENT SUBLIST OF GENES (50 one-genes);
  sum-genes                 ;COUNT THE SUM OF GENES MUTATIONS;
  apc-P-exponent            ;EXPONENT FOR THE APC PROBABILITY OF MUTATION 1 x 10 ^ - (gene-P-exponent);
  Kras-P-exponent           ; // ;
  p53-P-exponent            ; // ;
  other-mut-P-exponent      ; // ;
  born?                     ;WHEN A CELL IS CREATED SET ITS STATUS IN BORN? "TRUE";
  own-spotx
  own-spoty
  one-genes
  ID
  offspring
  ancestor
  descendant
  trait
  trait-d
  Z
  generation
  parent-who
  degree
]

TO SETUP

  random-seed 100

  ca

  ask patches with [pycor <= (growth-cell-radius * max-pycor / 100)] [set pcolor blue + 3]
  ask patches with [pycor > (growth-cell-radius * max-pycor / 100)]  [set pcolor green + 2]


  set n-repetitions 100

  let yy 1
  ;FILLS THE CRYPT WITH CELLS TRANSIT AND DIFFERENTIATED;
  repeat n-repetitions
  [
    let xx 1
    create-cells 25 ;CELLS' NUMBER PER LINE
    [
       set size cells-size
       set heading 0
       setxy xx yy
       set xx xx + 2
       if (pcolor = blue + 3)
       [
          set tipology "TRANSIT-AMPLIFYING"
          set color blue - 2 + random-float 3
          set mitosis-mean-time 24
          set mitosis-time -50
          while [mitosis-time < 16 or mitosis-time > max-cells-life-time]
          [set mitosis-time round random-normal mitosis-mean-time 2]
          ifelse (random 2 = 0) [set shape "realtransit"] [set shape "realtransit2"]
       ]
       if (pcolor = green + 2)
       [
          set tipology "DIFFERENTIATED"
          set mitosis-mean-time 30
          set color green - 2 + random-float 3 set mitosis-time -10
          ifelse (random 2 = 0) [set shape "realdiff"] [set shape "realdiff2"]
       ]
       set life-time (round (yy / 3) + ifelse-value (random 2 = 1)[random 10] [random -10])
       if (life-time < 0) [set life-time 0]
       set own-line yy
       set pre-adenoma? false
       set adenoma? false
       set tumoral? false
       set mitosis? false
       set B-block? false
       set p53-control false
       set born? false
       set APC [0 0]
       set Kras [0 0]
       set p53 [0 0]
       set one-genes [0 0]
       set genes n-values 5 [one-genes];sequenza di 50 geni
       set sum-genes 0
       set B-cat (1 - yy / 300)
       set apc-P-exponent apc-P
       set Kras-P-exponent Kras-P
       set p53-P-exponent p53-P
       set other-mut-P-exponent other-mut-P
       set offspring 0
       set ancestor false
       set descendant false
       set trait 0
      set parent-who random 2500
      if (parent-who = [parent-who] of one-of other cells)
      [
        set parent-who random 2500
      ]
      set degree 0

    ]
    set yy yy + 3
  ]
  set max-cells-line max [own-line] of cells

  set hours 0
  set n-cells-dead 0
  set n-hypoxia-death 0
  set met-death 0
  set energy-death 0


  set APC-mutations 0
  set Kras-mutations 0
  set p53-mutations 0
  set number-other-mut 0


;set-current-plot "LIFE-TIME-HISTOGRAM"
;set-plot-x-range
;precision (min [life-time] of cells) 1
;precision (max [life-time] of cells + 0.01) 1
;histogram [life-time] of cells
;
;set-current-plot "MITOSIS-TIME-HISTOGRAM"
;set-plot-x-range
;precision (min [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"]) 1
;precision (max [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"] + 0.01) 1
;set-plot-y-range 0 30
;set-current-plot-pen "TRANSIT"
;histogram [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"]

END

;STARTS SYSTEM DYNAMIC;
TO START

set hours hours + 1


  ask cells [

    set life-time life-time + 1 EFFECTS-OF-GENES-MUTATION

    if (mitosis? = true) [ set own-spotx pxcor set own-spoty pycor ]

    ifelse (ycor < 3) [set ID false] [set ID true]

  ]


;DEATH OR INACTIVATION OPERATIONS;
DEATH-INACTIVITY

set cells-line (max-cells-line - cells-size)

no-display
;REPEAT ALL THE PROCEDURES IN THE BOX;
repeat (n-repetitions - 1)
[

   ;1: NORMAL CELLS OF EVERY LINE MOVE FORWARD IF THEY FIND EMPTY SPACES;
   FORWARD-CELLS

   ;2: REPLICATION OPERATIONS FOR NORMAL CELLS WHO HAVE MOVED AND WITH MITOSIS ACTIVE;
   mitosis-MAIN-OPERATIONS

   ;3: MITOSIS FOR DIFFERENTIATED CELLS WHO FIND EMPTY SPACES;
   OTHER-mitosis-OPERATIONS

   ;4: FORCED MITOSIS FOR ALL CELL TYPES;
   FORCED-mitosis

   ;5: PRE-ADENOMA AND ADENOMA CELLS MOVE FORWARD MORE SLOWLY THAN THE OTHERS AND POSSIBLY GO DISORDERED MITOSIS;
   let adenoma-cells-line cells-line
   repeat 3
   [
     FORWARD-ADENOMA-CELLS adenoma-cells-line
     set adenoma-cells-line adenoma-cells-line + 1
   ]

   set cells-line cells-line - cells-size

]

;6: TUMORAL CELLS DON'T MOVE BUT GO IN DISORDERED MITOSIS;
mitosis-TUMORAL-CELLS

;7: PROCEDURE FOR hypoxia OF THE MICROENVIRONMENT AND THE CREATION OF IMMUNE SYSTEM CELLS;
ENVIRONMENT

let xx 1
let yy 1

;CREATE THE FIRST LINE OF CELLS FROM STAMINAL CELLS;

create-cells 25
[
   ;output-print xx
   set size cells-size
   set tipology "TRANSIT-AMPLIFYING"
   ifelse (random 2 = 0) [set shape "realtransit"] [set shape "realtransit2"]
   set color blue - 2 + random-float 3
   set life-time 0
   set mitosis-mean-time 24
   set mitosis-time -50
   while [mitosis-time < 16 or mitosis-time > max-cells-life-time]
   [set mitosis-time round random-normal mitosis-mean-time 2]
   set heading 0
   setxy xx yy
   set xx xx + 2
   set own-line yy
   set pre-adenoma? false
   set adenoma? false
   set tumoral? false
   set B-block? false
   set p53-control false
   set born? false
   ifelse (random-float 1 >= prob-stem-mut) [set APC [0 0]]
   [if-else (random 2 = 0) [set APC [1 0]] [set APC [0 1]]]

   ifelse (random-float 1 >= prob-stem-mut) [set Kras [0 0]]
   [if-else (random 2 = 0) [set KRAS [1 0]] [set KRAS [0 1]]]

   ifelse (random-float 1 >= prob-stem-mut) [set P53 [0 0]]
   [if-else (random 2 = 0) [set P53 [1 0]] [set P53 [0 1]]]

   set one-genes [0 0]
   set genes n-values 5 [one-genes]
   set sum-genes 0
   set B-cat (1 - yy / 300)
   set apc-P-exponent apc-P
   set Kras-P-exponent Kras-P
   set p53-P-exponent p53-P
   set other-mut-P-exponent other-mut-P
   set offspring 0
   set ancestor false
   set descendant false

    set trait 0

    set parent-who random 2500
      if (parent-who = [parent-who] of one-of other cells)
      [
        set parent-who random 2500
      ]
    set generation 0

]

;procedure di stop;
;if ((count cells with [tipology = "TRANSIT-AMPLIFYING"] / count cells) >= 0.75 or count cells = 10000 or

if (int(Hours / 24) = 365) [stop]

display

set count-cells (count cells with [ID = false])

;data-first-paper

;Price_eq_data

;DATA

;if (hours mod 120 = 0) [DATA_PRICE_EQ]

PLOTS

END

TO FORWARD-CELLS

   ask cells with [pre-adenoma? = false and own-line = cells-line]
   [
      if (not any? cells-on patch-ahead cells-size)
      [
         fd cells-size
         set own-line ycor
         set B-cat (1 - ycor / 300)

      REGULATORY-GENES

      if (tipology = "TRANSIT-AMPLIFYING" and B-block? = false and B-cat < random-float 0.80)
         [
            set tipology "DIFFERENTIATED"
            ifelse (random 2 = 0) [set shape "realdiff"] [set shape "realdiff2"]
            set color (green - 2) + random-float 3 ;set mitosis-time -10
         ]
         if (tipology = "TRANSIT-AMPLIFYING" and (life-time mod mitosis-time = 0))
         [set mitosis? true]
         if (tipology = "DIFFERENTIATED")
         [set mitosis? false]
      ]
   ]

END

TO FORWARD-ADENOMA-CELLS [acells-line]

   ask cells with [tumoral? = false and pre-adenoma? = true or adenoma? = true and own-line = acells-line]
   [

      ;if (not any? cells-on patch-right-and-ahead cells-size cells-size and not any? cells-on patch-left-and-ahead cells-size cells-size)

     ;[
      set degree random 90
      let target 0
      ifelse (random 2 = 0) [set target patch-right-and-ahead degree 1] [set target patch-left-and-ahead degree 1]

      move-to target
      set own-line ycor
      set B-cat (1 - ycor / 300)

      REGULATORY-GENES

      if (tipology = "TRANSIT-AMPLIFYING" and (life-time mod mitosis-time = 0))
      [
          set mitosis? true
          hatch-cells 1
          [
            set born? true
            set life-time 0
            move-to target set own-line ycor
            set offspring 0
          ]
          set offspring offspring + 1
          set mitosis? false
      ]


      if (tipology = "DIFFERENTIATED" and (life-time mod mitosis-time = 0)); and life-energy >= (life-energy / 3))
      [
         set mitosis? true
         hatch-cells 1
         [
           set born? true
           set life-time 0
           move-to target set own-line ycor
           set offspring 0

         ]
      set offspring offspring + 1
      set mitosis? false
      ]
    ]
   ;]


END


TO MITOSIS-MAIN-OPERATIONS

   ask cells with [tipology = "TRANSIT-AMPLIFYING" and tumoral? = false and own-line = (cells-line + cells-size) and mitosis? = true]
   [
      ; EVERY CELL CONTROL FOR EMPTY SPACES OR NON TUMORAL CELLS IN THE RIGHT AND IN THE LEFT;
      let right-cell 0
      if (any? cells-at 2 0)
      [
         set right-cell one-of cells-at 2 0
         if ([tumoral?] of right-cell = true or [sum p53] of right-cell = 2)[set right-cell 1]
      ]
      let left-cell 0
      if (any? cells-at -2 0)
      [
         set left-cell one-of cells-at -2 0
         if ([tumoral?] of left-cell = true or [sum p53] of left-cell = 2)[set left-cell 1]
      ]

      ;MULTIPLE IF-ELSE PROCEDURE;
      (ifelse
      right-cell = 0 or left-cell = 0

      ;IF THERE'S AT LEAST AN EMPTY SPACE IN THE RIGHT OR IN THE LEFT, THE CELL "ASKED" GOES IN MITOSIS;
      [
          GENES-MUTATION
          let xxx xcor
          hatch-cells 1
          [
             set born? true
             set life-time 0
             ifelse (right-cell = 0)
             [setxy (xxx + 2) own-line]
             [setxy (xxx - 2) own-line]
             set offspring 0

          ]
          set offspring offspring + 1
          set mitosis? false

      ]
      right-cell = 1 and left-cell = 1

      ;EITHER THE RIGHT AND THE LEFT CELLS ARE IMMORTAL, SO THE "ASKED" CELL CAN'T GOES IN MITOSIS;
      [
         set mitosis? false
      ]
      right-cell != 1 and left-cell = 1

      ;ONLY THE RIGHT CELL IS NOT IMMORTAL, SO THE "ASKED" CELL ELIMINATES IT AND GOES IN MITOSIS;
      [
         GENES-MUTATION
         ask right-cell [set n-cells-dead n-cells-dead + 1 die]
         let xxx xcor
         hatch-cells 1
         [
            set born? true
            set life-time 0
            setxy (xxx + 2) own-line
            set offspring 0

         ]
        set offspring offspring + 1
         set mitosis? false

       ]
       right-cell = 1 and left-cell != 1

       ;ONLY THE LEFT CELL IS NOT IMMORTAL, SO THE "ASKED" CELL ELIMINATES IT AND GOES MITOSIS
       [
          if (left-cell != 0)
          [
             GENES-MUTATION
             ask left-cell [set n-cells-dead n-cells-dead + 1 die]
             let xxx xcor
             hatch-cells 1
             [
                set born? true
                set life-time 0
                setxy (xxx - 2) own-line
                set offspring 0

             ]
             set offspring offspring + 1
             set mitosis? false

          ]
       ]
       right-cell != 1 and left-cell != 1

       ;BOTH ADJACENT CELLS ARE NOT IMMORTAL, SO THE "ASKED" CELL ELIMINATES THE OLDER ONE (OR THE ONE WITH LOWER FITNESS) AND GOES IN MITOSIS;
       [

          ifelse ([life-time] of right-cell > [life-time] of left-cell)
          [
             GENES-MUTATION
             ask right-cell [set n-cells-dead n-cells-dead + 1 die]
             let xxx xcor
             hatch-cells 1
             [
                set born? true
                set life-time 0
                setxy (xxx + 2) own-line
                set offspring 0

             ]
             set offspring offspring + 1
             set mitosis? false

          ]
          [
             GENES-MUTATION
             ask left-cell [set n-cells-dead n-cells-dead + 1 die]
             let xxx xcor
             hatch-cells 1
             [
                set born? true
                set life-time 0
                setxy (xxx - 2) own-line
                set offspring 0

             ]
             set offspring offspring + 1
             set mitosis? false
          ]
       ])
    ]

END

TO OTHER-MITOSIS-OPERATIONS

   ask patches with [pxcor mod 2 != 0 and pcolor = green + 2 and pycor = (cells-line + cells-size) and not any? cells-here]
   [
      ifelse (random 2 = 1)
      [
         if (any? (cells-on patch-at 2 0) with [pre-adenoma? = false])
         [
            ask one-of cells-on patch-at 2 0
            [
               GENES-MUTATION
               let xxx xcor
               set mitosis? true
               hatch-cells 1
               [
                  set born? true
                  set life-time 0
                  setxy (xxx - 2) own-line
                  set offspring 0

               ]
              set offspring offspring + 1
              set mitosis? false
            ]
         ]
      ]
      [
         if (any? (cells-on patch-at -2 0) with [pre-adenoma? = false])
         [
            ask one-of cells-on patch-at -2 0
            [
               GENES-MUTATION
               let xxx xcor
               set mitosis? true
               hatch-cells 1
               [
                  set born? true
                  set life-time 0
                  setxy (xxx + 2) own-line
                  set offspring 0

               ]
               set offspring offspring + 1
               set mitosis? false
            ]
         ]
      ]
   ]

END

TO FORCED-MITOSIS

   ask patches with [pxcor mod 2 != 0 and pycor = (cells-line + cells-size) and not any? cells-here]
   [
      let xxx pxcor
      if (any? cells with [own-line = (cells-line + cells-size) and pre-adenoma? = false])
      [
         ask one-of cells with [own-line = (cells-line + cells-size) and pre-adenoma? = false]
         [
            GENES-MUTATION
            set mitosis? true
            hatch-cells 1
            [
               set born? true
               set life-time 0
               setxy xxx own-line
               set offspring 0

            ]
            set offspring offspring + 1
            set mitosis? false
         ]
      ]
   ]

END

TO MITOSIS-TUMORAL-CELLS

   ask cells with [tumoral? = true and (life-time mod mitosis-time = 0)]
   [
     set degree random 360
     set mitosis? true
     hatch-cells 1
     [
      set born? true
      set life-time 0
      move-to patch-at degree 1
      set offspring 0

     ]
     set offspring offspring + 1
     set mitosis? false
   ]

END

TO GENES-MUTATION

;ALL THE PROCEDURES FOR CALCULATE PROBABILITY OF MUTATION, AND REPLACE THE ITEMS OF THE STRINGS WHO REPRESENT GENE ALLELES;

      let APC-prob-mutation-1 (10 ^ (- apc-P-exponent))
      let APC-prob-mutation-2 (10 ^ (- apc-P-exponent + 1))
      let Kras-prob-mutation (10 ^ (- Kras-P-exponent))
      let p53-prob-mutation-1 (10 ^ (- p53-P-exponent))
      let p53-prob-mutation-2 (10 ^ (- p53-P-exponent + 1))


      if (sum APC = 0) [if (random-float 1 < APC-prob-mutation-1) [set APC replace-item 0 APC 1]]
      if (sum APC = 1) [if (random-float 1 < APC-prob-mutation-2) [set APC replace-item 1 APC 1 set APC-mutations (APC-mutations + 1)]]
      if (sum Kras = 0) [if (random-float 1 < Kras-prob-mutation) [set Kras replace-item 0 Kras 1 set Kras-mutations (Kras-mutations + 1)]]
      if (sum Kras = 1) [if (random-float 1 < Kras-prob-mutation) [set Kras replace-item 1 Kras 1]]
      if (sum p53 = 0) [if (random-float 1 < p53-prob-mutation-1) [set p53 replace-item 0 p53 1]]
      if (sum p53 = 1) [if (random-float 1 < p53-prob-mutation-2) [set p53 replace-item 1 p53 1 set p53-mutations (p53-mutations + 1)]]



END

to REGULATORY-GENES

  let other-mut-Prob-mutation (10 ^ (- other-mut-P-exponent))


      foreach genes
      [
         ggenes ->

         if (sum ggenes = 0)
         [if (random-float 1 < other-mut-Prob-mutation)[set ggenes replace-item 0 ggenes 1 set sum-genes sum-genes + 1]]
         if (sum ggenes = 1)
         [if (random-float 1 < other-mut-Prob-mutation)[set ggenes replace-item 1 ggenes 1 set sum-genes sum-genes + 1]]

      ]

end

TO EFFECTS-OF-GENES-MUTATION

ifelse (tumoral? = false)

  [
;SUPPRESSION EFFECT OF APC;
      if (sum APC = 2)
      [
         set pre-adenoma? true
         set B-block? true
      ]

;SUPPRESSION EFFECT OF P53;
      if (sum p53 = 2)
      [
        if (APC-P-exponent > 1)
        [set APC-P-exponent APC-P-exponent - 0.5]
        if (Kras-P-exponent > 0)
        [set Kras-P-exponent Kras-P-exponent - 0.5]
        if (p53-P-exponent > 1)
        [set p53-P-exponent p53-P-exponent - 0.5]
        if (other-mut-P-exponent > 0)
        [set other-mut-P-exponent other-mut-P-exponent - 0.5]
      ]

;IPER-EXPRESSION OF K-RAS;
      if (sum Kras > 0)
      [
         if (mitosis-mean-time > 12) [set mitosis-mean-time 12]
         set mitosis-time -50
         while [mitosis-time < 5 or mitosis-time > max-cells-life-time]
         [set mitosis-time round random-normal mitosis-mean-time 2]
      ]

;HOW CELLS' FENOTYPE CHANGE WITH MUTATIONS;
      if (pre-adenoma? = true and adenoma? = false)
      [
        set color blue + 1
        set trait 1
      ]

      if (pre-adenoma? = true and sum Kras > 0)
      [
        set adenoma? true
        set color orange + sum-genes
        set trait 2
      ]

      if (adenoma? = true and sum p53 = 2)
      [
        set tumoral? true

         if (mitosis-mean-time > 10) [set mitosis-mean-time 10]
         set mitosis-time -50
         while [mitosis-time < 5 or mitosis-time > max-cells-life-time]
         [set mitosis-time round random-normal mitosis-mean-time 2]

        set shape "tumoral" set color red + sum-genes
        set trait 3
      ]
]

;SUM-GENES REACH MAXIMUM 100;
  [
    if (tumoral? = true and sum-genes = 100)
    [set shape "tumoral2" set color red + sum-genes]
  ]

;DEREGULATION DUE TO SUM-GENES

    if (sum-genes > 25)
    [
      if (APC-P-exponent > 1)
        [set APC-P-exponent APC-P-exponent - 0.5]
        if (Kras-P-exponent > 0)
        [set Kras-P-exponent Kras-P-exponent - 0.5]
        if (p53-P-exponent > 1)
        [set p53-P-exponent p53-P-exponent - 0.5]
        if (other-mut-P-exponent > 0)
        [set other-mut-P-exponent other-mut-P-exponent - 0.5]
    ]



END

TO ENVIRONMENT

  ;PATCHES WITH hypoxiaITY OVER THE THRESHOLD PRODUCE A VARIABLE NUMBER OF IMMUNE-SYSTEM CELLS;
    ask patches
  [
    set hypoxia (count cells-here with [pre-adenoma? = true or adenoma? = true or tumoral? = true])
    if (hypoxia >= immune-threshold)
    [
       sprout-immune-cells immune-system
       [
         set shape "lymph"
         set size cells-size
         set task count n-of 1 cells
       ]
    ]
  ]

  ;IF IMMUNE-SYSTEM CELLS DO NOT FIND ANY DANGEROUS CELL IN RADIUS 1 THEY DIE, OTHERWISE THEY KILL THEM;
  ask immune-cells
  [
    fd 2
    set immune-life immune-life + 1
    if (not any? cells in-radius 2 with [pre-adenoma? = true or adenoma? = true or tumoral? = true])
    [set heading random 360 fd 1]
    if (any? cells in-radius 2 with [pre-adenoma? = true or adenoma? = true or tumoral? = true])
    [
        let one-cell one-of cells with [pre-adenoma? = true or adenoma? = true or tumoral? = true]
        move-to one-cell
        ask one-cell [die]
        set task task + 1
        set immune-death immune-death + 1
    ]
  ]



END


TO DEATH-INACTIVITY

;DEATH CLAUSES;
;a)
  ask cells with [tumoral? = false and (life-time > max-cells-life-time)] [set born? false set n-cells-dead n-cells-dead + 1 die]
;b)
  ask cells with [tumoral? = false and ycor > (max-cells-line - cells-size)] [set born? false die]
;c)
  ask cells with [tumoral? = true and (life-time > max-tumoral-cells-life-time)][set born? false set n-tumoral-dead n-tumoral-dead + 1 die]
;d)
  ask cells [if [hypoxia] of patch-here >= hypoxia-threshold [set born? false set n-hypoxia-death n-hypoxia-death + 1 die]]
;e)
  ask cells [if (sum-genes = 5 and sum p53 < 2) [set born? false set p53-control true set p53-dead p53-dead + 1 die]]
;f)
  ask immune-cells [if (task >= number-task) or immune-life > 10 [die]]
;g)
  ask cells with [pre-adenoma? = true or adenoma? = true or tumoral? = true and ycor > (max-cells-line - cells-size)] [set born? false die]

END

TO PLOTS
display
;PLOTS FOR CONTROL THE RELEVANT PARAMETRES OF THE SIMULATION;

set-current-plot "CELLS-GROWTH-PLOT"
set-current-plot-pen "TRANSIT"
plotxy hours count cells with [tipology = "TRANSIT-AMPLIFYING" and pre-adenoma? = false and ID = true]
set-current-plot-pen "DIFFERENT"
plotxy hours count cells with [tipology = "DIFFERENTIATED" and pre-adenoma? = false and ID = true]
set-current-plot-pen "PRE-ADENOMA"
plotxy hours count cells with [pre-adenoma? = true and adenoma? = false and ID = true]
set-current-plot-pen "ADENOMA"
plotxy hours count cells with [adenoma? = true and tumoral? = false and ID = true]
set-current-plot-pen "TUMORAL"
plotxy hours count cells with [tumoral? = true and ID = true]
set-current-plot-pen "IMMUNE"
plotxy hours count immune-cells

;set-current-plot "Pre-prey"
;set-current-plot-pen "1"
;  histogram [parent-who] of cells with [ID = true]
;plotxy count immune-cells count cells with [pre-adenoma? = true and ID = true]
;
;  set-current-plot "WHO"
;  set-current-plot-pen "w"
;  histogram [parent-who] of cells with [ID = true]
;
;set-current-plot "CELLS-DEATH-PLOT"
;set-current-plot-pen "n-death"
;plotxy hours round (n-cells-dead / hours)
;set-current-plot-pen "n-hypoxia"
;plotxy hours round (n-hypoxia-death / hours)
;set-current-plot-pen "immune-death"
;plotxy hours round (immune-death / hours)
;set-current-plot-pen "tum-death"
;plotxy hours round (n-tumoral-dead / hours)
; set-current-plot-pen "sole"
;plotxy hours round (cells-dead / hours)
;
;set-current-plot "LIFE-TIME-HISTOGRAM"
;set-plot-x-range
;precision (min [life-time] of cells) 1
;precision (max [life-time] of cells + 0.01) 1
;histogram ([life-time] of cells)
;
;set-current-plot "mitosis-TIME-HISTOGRAM"
;set-plot-x-range
;precision (min [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"]) 1
;precision (max [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"] + 0.01) 1
;set-plot-y-range 0 30
;set-current-plot-pen "TRANSIT"
;histogram [mitosis-time] of cells with [tipology = "TRANSIT-AMPLIFYING"]
;set-current-plot-pen "DIFF"
;histogram [mitosis-time] of cells with [tipology = "DIFFERENTIATED"]
;
set-current-plot "HYPOXIA"
set-plot-y-range 0 100
set-current-plot-pen "HYP"
histogram [hypoxia] of patches
;
;set-current-plot "proliferation"
;set-plot-x-range
;(min [offspring] of cells) (max [offspring] of cells)
;set-plot-y-range
;0 count cells with [offspring = 0]
;set-current-plot-pen "normal"
;histogram [ offspring ] of cells with [pre-adenoma? = false]
;set-current-plot-pen "tumoral"
;histogram [ offspring ] of cells with [tumoral? = true]
;set-current-plot-pen "adenoma"
;histogram [ offspring ] of cells with [adenoma? = true and tumoral? = false]
;set-current-plot-pen "pre-ade"
;histogram [ offspring ] of cells with [pre-adenoma? = true and adenoma? = false]



;set-current-plot "genotype-frequency"
;set-plot-y-range 0 1.5
;
;  set-current-plot-pen "apc0"
;plot ((count cells with [sum APC = 0]) / (count cells))
;;------------------------------APC heterozygotes------------------------------
;  set-current-plot-pen "apc1-"
;plot ((count cells with [item 0 APC = 1 and item 1 APC = 0]) / (count cells))
;  set-current-plot-pen "apc-1"
;plot ((count cells with [item 1 APC = 1 and item 0 APC = 0]) / (count cells))
;;-----------------------------------------------------------------------------
;  set-current-plot-pen "apc2"
;plot ((count cells with [sum APC > 1]) / (count cells))
;
;  set-current-plot-pen "kras0"
;plot ((count cells with [sum kras = 0]) / (count cells))
;;------------------------------Kras heterozygotes------------------------------
;  set-current-plot-pen "kras1-"
;plot ((count cells with [item 0 kras = 1 and item 1 kras = 0]) / (count cells))
;  set-current-plot-pen "kras-1"
;plot ((count cells with [item 1 kras = 1 and item 0 kras = 0]) / (count cells))
;;------------------------------------------------------------------------------
;  set-current-plot-pen "kras2"
;plot ((count cells with [sum kras > 1]) / (count cells))
;
;  set-current-plot-pen "p530"
;plot ((count cells with [sum p53 = 0]) / (count cells))
;;------------------------------P53 heterozygotes-------------------------------
;  set-current-plot-pen "p531-"
;plot ((count cells with [item 0 p53 = 1 and item 1 p53 = 0]) / (count cells))
;
;   set-current-plot-pen "p53-1"
;plot ((count cells with [item 1 p53 = 1 and item 0 p53 = 0]) / (count cells))
;;------------------------------------------------------------------------------
;  set-current-plot-pen "p532"
;plot ((count cells with [sum p53 > 1]) / (count cells))


END

to data-first-paper


;  file-open "cell-patch.csv"
;
;  file-print ""
;
;  file-write [hypoxia] of patches file-print ""
;
;  file-close


;  file-open "population.csv"
;
;  file-print ""
;
;  file-write count cells with [tipology = "TRANSIT-AMPLIFYING" and pre-adenoma? = false and ID = true]
;
;  file-write count cells with [tipology = "DIFFERENTIATED" and pre-adenoma? = false and ID = true]
;
;  file-write count cells with [pre-adenoma? = true and adenoma? = false and ID = true]
;
;  file-write count cells with [adenoma? = true and tumoral? = false and ID = true]
;
;  file-write count cells with [tumoral? = true and ID = true]
;
;  file-write count immune-cells
;
;  file-close


;  file-open "Death.csv"
;
;  file-print ""
;
;  file-write count cells with [tipology = "TRANSIT-AMPLIFYING" and pre-adenoma? = false and ID = true]
;
;  file-write count cells with [tipology = "DIFFERENTIATED" and pre-adenoma? = false and ID = true]
;
;  file-write count cells with [pre-adenoma? = true and adenoma? = false and ID = true]
;
;  file-write count cells with [adenoma? = true and tumoral? = false and ID = true]
;
;  file-write count cells with [tumoral? = true and ID = true]
;
;  file-write count immune-cells
;
;  file-close


;  file-open "genotype.csv"
;
;  file-print ""
;
;  file-write ((count cells with [sum APC = 0]) / (count cells))
;  file-write ((count cells with [item 0 APC = 1 and item 1 APC = 0]) / (count cells))
;  file-write ((count cells with [item 0 APC = 0 and item 1 APC = 1]) / (count cells))
;  file-write ((count cells with [sum APC > 1]) / (count cells))
;  ;----------------------------------------------------------------------------------
;  file-write ((count cells with [sum kras = 0]) / (count cells))
;  file-write ((count cells with [item 0 kras = 1 and item 1 kras = 0]) / (count cells))
;  file-write ((count cells with [item 0 kras = 0 and item 1 kras = 1]) / (count cells))
;  file-write ((count cells with [sum kras > 1]) / (count cells))
;  ;------------------------------------------------------------------------------------
;  file-write ((count cells with [sum p53 = 0]) / (count cells))
;  file-write ((count cells with [item 0 p53 = 1 and item 1 p53 = 0]) / (count cells))
;  file-write ((count cells with [item 1 p53 = 1 and item 0 p53 = 0]) / (count cells))
;  file-write ((count cells with [sum p53 > 1]) / (count cells))
;  ;------------------------------------------------------------------------------------
;
;  ask cells with [ID = true]
;  [
;
;     file-write [offspring] of self file-write [parent-who] of self
;
;       file-write (count immune-cells / count cells) file-write (n-hypoxia-death / count cells)
;
;    file-print ""
;  ]
;
;  file-close




;  file-open "P53 death.csv"
;
;  file-print ""
;
;  file-write round (p53-dead / hours) file-close
;
;
;  file-open "sum-genes.csv"
;
;  file-print ""
;
;  ask cells with [ID = true]  [file-write [sum-genes] of self file-print ""]
;
;  file-close



end



to Price_eq_data

  file-open "Price_eq_data.csv"

  file-print ""

  ask cells with [ID = true] [

    file-write sum [apc] of self file-write sum [kras] of self file-write sum [p53] of self

     file-write sum [apc] of self + sum [kras] of self + sum [p53] of self

      file-write [offspring] of self file-write [parent-who] of self

       file-write (count immune-cells / count cells) file-write (n-hypoxia-death / count cells)

    file-print ""]

;  file-close
;
;  file-open "environment.csv"
;
;  file-print ""
;
;  file-write count immune-cells file-write round (n-hypoxia-death / hours)

  file-close

end


to DATA_PRICE_EQ

;  file-open "ancestor sample.csv"
;
;  file-print ""
;
;  ask cells with [ID = true] [file-write Z file-write trait file-write [offspring] of self file-write hours file-print ""]
;
;  file-close

  file-open "descendant sample.csv"

  file-print ""

  ask cells with [ID = true and descendant = true] [file-write Z file-write trait file-write [offspring] of self file-write hours file-print ""]

  file-close

end

to DATA


;  file-open "Cell death.csv"
;
;  file-print ""
;
;  file-write round (n-cells-dead / hours)
;
;  file-write round (n-hypoxia-death / hours)
;
;  file-write round (immune-death / hours)
;
;  file-write round (n-tumoral-dead / hours)
;
;  file-close
;
;
;
;  file-open "Immune system population.csv"
;
;  file-print ""
;
;  file-write (count immune-cells / hours)
;
;  file-write (count cells with [ID = true and pre-adenoma? = true and adenoma? = false] / hours)
;
;  file-write (count cells with [ID = true and pre-adenoma? = true and adenoma? = true and tumoral? = false] / hours)
;
;  file-write (count cells with [ID = true and pre-adenoma? = true and adenoma? = true and tumoral? = true] / hours)
;
;  file-close
;
;
;
;
;
;
;  file-open "Price-data-gene.csv"
;
;  file-print ""
;  ;------------------frequencies of normal cells---------------------------------------------
;
;  file-write  ((count cells with [sum APC = 0 and sum kras = 0 and sum p53 = 0]) / (count cells))
;
;  ;------------------frequencies of pre-adenoma cells---------------------------------------------
;
;  file-write  ((count cells with [sum APC = 2 and sum kras = 0 and sum p53 = 0]) / (count cells))
;
;  ;------------------frequencies of adenoma cells---------------------------------------------
;
;  file-write  ((count cells with [sum APC = 2 and sum kras > 0 and sum p53 = 0]) / (count cells))
;
;  ;------------------frequencies of tumoral cells---------------------------------------------
;
;  file-write  ((count cells with [sum APC = 2 and sum kras > 0 and sum p53 = 2]) / (count cells))
;
;  file-close
;
;
;  file-open "Price-data-fitness.csv"
;
;  file-print ""
;
;  ;----------------fitness normal cells------------------------------------------------------------
;
;  file-write  (sum [offspring] of cells with [sum APC = 0 and sum kras = 0 and sum p53 = 0]) / (count cells with [ID = true and sum APC = 0 and sum kras = 0 and sum p53 = 0] + 0.0001)
;
;  ;----------------fitness pre-adenoma cells------------------------------------------------------------
;
;  file-write  (sum [offspring] of cells with [sum APC = 2 and sum kras = 0 and sum p53 = 0]) / (count cells with [ID = true and sum APC = 2 and sum kras = 0 and sum p53 = 0] + 0.0001)
;
;  ;----------------fitness adenoma cells------------------------------------------------------------
;
;  file-write  (sum [offspring] of cells with [sum APC = 2 and sum kras > 0 and sum p53 = 0]) / (count cells with [ID = true and sum APC = 2 and sum kras > 0 and sum p53 = 0] + 0.0001)
;
;  ;----------------fitness tumoral cells------------------------------------------------------------
;
;  file-write  (sum [offspring] of cells with [sum APC = 2 and sum kras > 0 and sum p53 = 2]) / (count cells with [ID = true and sum APC = 2 and sum kras > 0 and sum p53 = 2] + 0.0001)
;
;
;  file-close
;

;  file-open "Allele-data.csv"
;
;  file-print ""
;  ;------------------frequencies of APC alleles---------------------------------------------
;  file-write  ((count cells with [sum APC = 0]) / (count cells))
;
;  file-write ((count cells with [sum APC = 1]) / (count cells))
;
;  file-write ((count cells with [sum APC = 2]) / (count cells))
;
;  ;------------------frequencies of KRAS alleles---------------------------------------------
;  file-write  ((count cells with [sum Kras = 0]) / (count cells))
;
;  file-write ((count cells with [sum kras = 1]) / (count cells))
;
;  file-write ((count cells with [sum Kras = 2]) / (count cells))
;
;  ;------------------frequencies of P53 alleles---------------------------------------------
;  file-write  ((count cells with [sum p53 = 0]) / (count cells))
;
;  file-write ((count cells with [sum p53 = 1]) / (count cells))
;
;  file-write ((count cells with [sum p53 = 2]) / (count cells))
;
;  ;-----------------------------------------------------------------------------------------
;
;  file-write (count immune-cells / count cells)
;
;  file-write (n-hypoxia-death / count cells)
;
;  ;-----------------------------------------------------------------------------------------
;
;  file-close
;
;
;;-------------------------------------------------------------------------------------------------
;
; file-open "Fitness.csv"
;
;  ;------APC-----------------
;  file-print ""
;
;  file-write (sum [offspring] of cells with [ID = true and sum APC = 0]) / (count cells with [ID = true and sum APC = 0] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum APC = 1]) / (count cells with [ID = true and sum APC = 1] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum APC = 2]) / (count cells with [ID = true and sum APC = 2] + 0.0001)
;
;  ;------kras----------------
;
;  file-write (sum [offspring] of cells with [ID = true and sum kras = 0]) / (count cells with [ID = true and sum kras = 0] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum kras = 1]) / (count cells with [ID = true and sum kras = 1] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum kras = 2]) / (count cells with [ID = true and sum kras = 2] + 0.0001)
;
;  ;------p53-----------------
;
;   file-write (sum [offspring] of cells with [ID = true and sum p53 = 0]) / (count cells with [ID = true and sum p53 = 0] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum p53 = 1]) / (count cells with [ID = true and sum p53 = 1] + 0.0001)
;
;  file-write (sum [offspring] of cells with [ID = true and sum p53 = 2]) / (count cells with [ID = true and sum p53 = 2] + 0.0001)
;
;
;  file-close



  file-open "Prol-ycor.csv" ask cells with [mitosis? = true] [ file-write own-spotx  file-write own-spoty file-write offspring file-write Z file-print ""]
  file-close

end
@#$#@#$#@
GRAPHICS-WINDOW
278
10
432
895
-1
-1
2.92
1
10
1
1
1
0
1
0
1
0
49
0
299
0
0
1
ticks
30.0

BUTTON
6
10
197
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
1
430
179
463
growth-cell-radius
growth-cell-radius
0
96
36.0
1
1
NIL
HORIZONTAL

SLIDER
1
466
179
499
cells-size
cells-size
0
30
3.0
1
1
NIL
HORIZONTAL

SLIDER
0
359
179
392
max-cells-life-time
max-cells-life-time
0
100
100.0
1
1
NIL
HORIZONTAL

BUTTON
198
10
279
141
NIL
start
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
438
10
869
232
CELLS-GROWTH-PLOT
hours
cells
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"TRANSIT" 1.0 0 -13791810 true "" ""
"DIFFERENT" 1.0 0 -14835848 true "" ""
"PRE-ADENOMA" 1.0 0 -4079321 true "" ""
"ADENOMA" 1.0 0 -955883 true "" ""
"TUMORAL" 1.0 0 -2674135 true "" ""
"IMMUNE" 1.0 0 -2064490 true "" ""

PLOT
714
395
965
563
LIFE-TIME-HISTOGRAM
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

MONITOR
90
95
195
140
NIL
hours
5
1
11

MONITOR
92
240
173
285
mean-lifetime
mean [life-time] of cells with [tipology != \"STEM\"]
1
1
11

MONITOR
5
48
87
93
# cells
count cells with [ID = true]
0
1
11

SLIDER
1
913
135
946
waiting-time
waiting-time
0
1
0.0
0.01
1
NIL
HORIZONTAL

PLOT
438
235
869
393
CELLS-DEATH-PLOT
hours
n-cells/cells
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"n-death" 1.0 0 -16777216 true "" ""
"n-hypoxia" 1.0 0 -13840069 true "" ""
"immune-death" 1.0 0 -2064490 true "" ""
"tum-death" 1.0 0 -5825686 true "" ""
"sole" 1.0 0 -6459832 true "" ""

MONITOR
6
95
86
140
NIL
n-cells-dead
0
1
11

PLOT
438
395
711
562
MITOSIS-TIME-HISTOGRAM
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"TRANSIT" 1.0 0 -13345367 true "" ""
"DIFF" 1.0 0 -14439633 true "" ""

MONITOR
90
48
195
93
days
int (hours / 24)
0
1
11

SLIDER
0
609
178
642
apc-P
apc-P
1
9
4.0
1
1
NIL
HORIZONTAL

SLIDER
0
643
178
676
Kras-P
Kras-P
1
9
4.0
1
1
NIL
HORIZONTAL

SLIDER
0
678
178
711
p53-P
p53-P
1
9
4.0
1
1
NIL
HORIZONTAL

SLIDER
0
876
135
909
methilation
methilation
0
1
0.0
0.1
1
NIL
HORIZONTAL

SLIDER
0
395
179
428
max-tumoral-cells-life-time
max-tumoral-cells-life-time
0
150
150.0
1
1
NIL
HORIZONTAL

SLIDER
0
713
178
746
other-mut-P
other-mut-P
1
9
4.0
1
1
NIL
HORIZONTAL

MONITOR
6
190
88
235
preadenoma
count cells with [pre-adenoma?]
0
1
11

MONITOR
6
239
88
284
adenoma
count cells with [adenoma?]
0
1
11

MONITOR
5
288
86
333
tumoral
count cells with [tumoral?]
0
1
11

SLIDER
0
502
179
535
hypoxia-threshold
hypoxia-threshold
0
5000
10.0
1
1
NIL
HORIZONTAL

SLIDER
0
539
178
572
immune-system
immune-system
0
10
5.0
1
1
NIL
HORIZONTAL

PLOT
968
395
1252
563
HYPOXIA
NIL
NIL
1.0
10.0
0.0
10.0
true
false
"" ""
PENS
"HYP" 1.0 1 -8431303 true "" ""

PLOT
872
11
1363
394
GENOTYPE-FREQUENCY
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"apc0" 1.0 0 -13345367 true "" ""
"apc2" 1.0 2 -13791810 true "" ""
"kras1-" 1.0 0 -955883 true "" ""
"kras2" 1.0 2 -6459832 true "" ""
"p530" 1.0 0 -15302303 true "" ""
"p532" 1.0 2 -13840069 true "" ""
"apc1-" 1.0 0 -11221820 true "" ""
"p531-" 1.0 0 -10899396 true "" ""
"kras0" 1.0 0 -4699768 true "" ""
"apc-1" 1.0 2 -8990512 true "" ""
"kras-1" 1.0 2 -955883 true "" ""
"p53-1" 1.0 2 -12087248 true "" ""

SLIDER
0
573
178
606
number-task
number-task
0
50
10.0
1
1
NIL
HORIZONTAL

PLOT
438
564
966
804
PROLIFERATION
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"normal" 1.0 0 -16777216 true "" ""
"tumoral" 1.0 0 -2674135 true "" ""
"adenoma" 1.0 0 -6459832 true "" ""
"pre-ade" 1.0 0 -13345367 true "" ""

MONITOR
90
289
195
334
mitosis-mean-time
mean [mitosis-mean-time] of cells
1
1
11

SLIDER
0
749
178
782
prob-stem-mut
prob-stem-mut
0
1
0.5
0.1
1
NIL
HORIZONTAL

MONITOR
92
191
149
236
trans
count cells with [tipology = \"TRANSIT-AMPLIFYING\"] / count cells
2
1
11

MONITOR
145
191
202
236
diff
count cells with [tipology = \"DIFFERENTIATED\"] / count cells
2
1
11

PLOT
968
569
1281
795
Pre-prey
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"1" 1.0 1 -16777216 true "" ""

PLOT
1256
399
1416
519
WHO
NIL
NIL
0.0
1000.0
0.0
10.0
true
false
"" ""
PENS
"W" 1.0 1 -2674135 true "" ""

SLIDER
0
785
179
818
immune-threshold
immune-threshold
0
100
5.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## COLORECTAL CANCER DEVELOPMENT (CRC)

CRC is an excellent system to simulate for two reasons in particular: first of all, it is a type of tumor of which we know a part of the multi-step progression of genetic mutations that favor the tumor phenotype; secondly, if we consider the cells immerged in their microenvironment as individuals of a population that changes over time, we have a model in which observe evolutionary mechanisms 
(speciation, selection of individuals, adaptation, etc.) 
on temporal scales accessible to the human investigation.

## HOW IT WORKS

The colonic crypt, modeled with the ABM software Netlogo, is displayed in a 2D space filled with two types of cells: transit-amplifying and differentiated. These cells move upward towards the crypt mouth where they come out the simulation, and, during the movement they replicate following given parameters. The cell’s behavior can be handled allowing experimentation of various scenarios (e.g., probability of gene mutation,
maximum lifespan, mitosis time, etc.) and it’s possible to follow the simulation progress using plots and monitors. The simulation starts with physiological features that change over time in function of genetic driver mutations; once that mutations occur, environmental and populational factors influence the cells behavior too.

## HOW TO USE IT

The model interface presents several sliders, some control driver genes probability of mutations (APC, K-ras, P53) and a group of genes with no specific functions, except as threshold for genome instability and for subclonal differentiation; others control maximum lifetime, the acidity of the micro-environment, the amount of immune system cells. 
Aside that, it’s possible to change parameters directly through the code thus increasing
the model’s potential. Given that, one can try the effect of increase the mutation rate of a single or more gene at a time or modify the number of immune cells produced to fight the cancer cells. In a more precise way, it’s possible to modify the cells movement, their mitosis time, causes of death, etc.

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

banana
false
0
Polygon -7500403 false true 25 78 29 86 30 95 27 103 17 122 12 151 18 181 39 211 61 234 96 247 155 259 203 257 243 245 275 229 288 205 284 192 260 188 249 187 214 187 188 188 181 189 144 189 122 183 107 175 89 158 69 126 56 95 50 83 38 68
Polygon -7500403 true true 39 69 26 77 30 88 29 103 17 124 12 152 18 179 34 205 60 233 99 249 155 260 196 259 237 248 272 230 289 205 284 194 264 190 244 188 221 188 185 191 170 191 145 190 123 186 108 178 87 157 68 126 59 103 52 88
Line -16777216 false 54 169 81 195
Line -16777216 false 75 193 82 199
Line -16777216 false 99 211 118 217
Line -16777216 false 241 211 254 210
Line -16777216 false 261 224 276 214
Polygon -16777216 true false 283 196 273 204 287 208
Polygon -16777216 true false 36 114 34 129 40 136
Polygon -16777216 true false 46 146 53 161 53 152
Line -16777216 false 65 132 82 162
Line -16777216 false 156 250 199 250
Polygon -16777216 true false 26 77 30 90 50 85 39 69

blood vessel
true
0
Polygon -2674135 true false 150 15 165 30 165 30 165 60 180 75 180 90 165 105 150 120 150 135 165 150 165 165 180 180 180 210 180 225 180 240 165 255 180 270 180 285 150 255 120 255 105 240 150 240 150 240 165 225 150 210 150 195 135 180 135 165 135 165 120 135 120 105 150 90 105 75 75 90 45 75 45 45 60 75 60 75 90 60 120 60 150 75 150 75 150 45 135 30 135 0

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

cell
true
10
Polygon -13345367 true true 150 15 105 15 75 45 75 255 105 285 150 285 195 285 225 255 225 45 225 45 195 15 150 15
Polygon -16777216 false false 105 15 75 45 75 255 105 285 195 285 225 255 225 45 195 15 105 15
Circle -11221820 true false 117 117 66
Circle -16777216 false false 105 15 0

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -16777216 true false 0 0 300
Circle -7500403 true true 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

enterocita
true
0
Polygon -13345367 true false 150 15 105 15 75 45 75 255 105 285 150 285 195 285 225 255 225 45 225 45 195 15 150 15
Polygon -16777216 false false 105 15 75 45 75 255 105 285 195 285 225 255 225 45 195 15 105 15
Circle -11221820 true false 117 117 66
Circle -16777216 false false 105 15 0

enterocita differenziato
true
0
Polygon -14835848 true false 150 15 105 15 75 45 75 255 105 285 150 285 195 285 225 255 225 45 225 45 195 15 150 15
Polygon -16777216 false false 105 15 75 45 75 255 105 285 195 285 225 255 225 45 195 15
Circle -11221820 true false 120 120 60
Circle -16777216 false false 105 15 0

enterocita mut
true
0
Polygon -2674135 true false 150 15 105 15 75 45 75 255 105 285 150 285 195 285 225 255 225 45 225 45 195 15 150 15
Polygon -16777216 false false 105 15 75 45 75 255 105 285 195 285 225 255 225 45 195 15
Circle -11221820 true false 120 120 60
Circle -16777216 false false 105 15 0

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

immune
false
0
Polygon -5825686 true false 150 15 15 120 60 285 240 285 285 120
Circle -1 true false 135 150 30
Circle -5825686 true false 105 0 90
Circle -5825686 true false 0 75 90
Circle -5825686 true false 30 210 90
Circle -5825686 true false 180 210 90
Circle -5825686 true false 210 75 90
Circle -1 true false 105 180 30
Circle -1 true false 165 180 30
Circle -1 true false 180 120 30
Circle -1 true false 135 120 30
Circle -1 true false 105 135 30
Circle -1 true false 105 90 30

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

lymph
false
5
Polygon -7500403 true false 135 30 105 30 90 45 60 75 45 105 45 135 45 180 60 225 75 255 120 270 195 270 225 240 240 225 255 195 255 180 255 165 240 90 195 45
Circle -2064490 true false 60 60 178

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

realdiff
false
5
Polygon -10899396 true true 165 30 120 30 75 60 45 90 45 135 45 150 45 195 60 225 75 255 105 285 180 285 210 270 225 240 240 210 255 180 255 135 240 60 195 30
Polygon -16777216 false false 120 30 75 60 45 90 45 195 75 255 105 285 180 285 210 270 255 180 255 135 240 60 195 30 120 30
Circle -1 true false 116 131 67

realdiff2
false
5
Polygon -10899396 true true 150 15 75 30 45 60 30 90 30 135 30 165 30 195 45 225 60 255 90 285 135 300 195 285 225 270 240 240 240 195 240 135 225 75 210 45
Circle -1 true false 101 116 67
Polygon -16777216 false false 75 30 45 60 45 60 30 90 30 195 60 255 90 285 135 300 195 285 225 270 240 240 240 135 225 75 225 75 210 45 150 15

realtransit
false
10
Polygon -13345367 true true 150 15 75 30 45 60 30 90 30 135 30 165 30 195 45 225 60 255 90 285 135 300 195 285 225 270 240 240 240 195 240 135 225 75 210 45
Polygon -16777216 false false 75 30 45 60 45 60 30 90 30 195 60 255 90 285 135 300 195 285 225 270 240 240 240 135 225 75 225 75 210 45 150 15
Circle -1 true false 101 116 67

realtransit2
false
10
Polygon -13345367 true true 165 30 120 30 75 60 45 90 45 135 45 150 45 195 60 225 75 255 105 285 180 285 210 270 225 240 240 210 255 180 255 135 240 60 195 30
Polygon -16777216 false false 120 30 75 60 45 90 45 195 75 255 105 285 180 285 210 270 255 180 255 135 240 60 195 30 120 30
Circle -1 true false 116 131 67

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

stem
true
0
Polygon -6459832 true false 150 15 105 15 75 45 75 255 105 285 150 285 195 285 225 255 225 45 225 45 195 15 150 15
Polygon -16777216 false false 105 15 75 45 75 255 105 285 195 285 225 255 225 45 195 15
Circle -11221820 true false 120 120 60
Circle -16777216 false false 105 15 0

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

tumoral
false
2
Polygon -955883 true true 150 30 150 30 90 45 45 75 45 165 60 210 60 270 105 285 180 285 225 255 255 180 225 90 225 45 210 30
Circle -1 true false 101 116 67
Polygon -16777216 false false 45 75 90 45 150 30 210 30 210 30 225 45 225 90 255 180 225 255 180 285 105 285 60 270 60 210 45 165

tumoral2
false
2
Polygon -955883 true true 150 30 120 45 90 75 90 90 75 135 75 150 75 210 105 270 135 285 165 285 195 255 210 195 210 150 210 90 195 60 180 45 150 30
Circle -1 true false 116 131 67
Polygon -16777216 false false 150 30 120 45 90 75 90 90 75 135 75 210 105 270 135 285 165 285 195 255 210 195 210 90 195 60 180 45

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="RUN-1" repetitions="10" runMetricsEveryStep="true">
    <setup>SETUP</setup>
    <go>START</go>
    <timeLimit steps="480"/>
    <metric>hours</metric>
    <metric>plot genotype-frequency</metric>
    <enumeratedValueSet variable="apc-P">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune-system">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kras-P">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-stem-mut">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-task">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="other-mut-P">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p53-P">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hypoxia-threshold">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
