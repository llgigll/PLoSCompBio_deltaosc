(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 9.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1063,         20]
NotebookDataLength[     18269,        625]
NotebookOptionsPosition[     17895,        587]
NotebookOutlinePosition[     18376,        608]
CellTagsIndexPosition[     18333,        605]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "1", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "Aie"], " ", 
     SubscriptBox["J", "ei"]}], "+", 
    RowBox[{
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["J", "Nie"]}]}], ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["J", "ii"]}], "+", 
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["J", "Nee"]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "2", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["J", "Nie"], " ", 
     RowBox[{"(", 
      RowBox[{
       SubscriptBox["\[Tau]", "Aee"], "+", 
       SubscriptBox["\[Tau]", "Aie"], "+", 
       SubscriptBox["\[Tau]", "ii"], "+", 
       SubscriptBox["\[Tau]", "Nee"]}], ")"}]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aie"], " ", 
     SubscriptBox["J", "ei"], " ", 
     RowBox[{"(", 
      RowBox[{
       SubscriptBox["\[Tau]", "Aee"], "+", 
       SubscriptBox["\[Tau]", "ii"], "+", 
       SubscriptBox["\[Tau]", "Nee"], "+", 
       SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["J", "Nee"], " ", 
     RowBox[{"(", 
      RowBox[{
       SubscriptBox["\[Tau]", "Aee"], "+", 
       SubscriptBox["\[Tau]", "Aie"], "+", 
       SubscriptBox["\[Tau]", "ei"], "+", 
       SubscriptBox["\[Tau]", "Nie"]}], ")"}]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["J", "ii"], " ", 
     RowBox[{"(", 
      RowBox[{
       SubscriptBox["\[Tau]", "Aie"], "+", 
       SubscriptBox["\[Tau]", "ei"], "+", 
       SubscriptBox["\[Tau]", "Nee"], "+", 
       SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "3", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["J", "Nie"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "ii"], " ", 
        SubscriptBox["\[Tau]", "Nee"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "ii"], "+", 
          SubscriptBox["\[Tau]", "Nee"]}], ")"}]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "Aie"], "+", 
          SubscriptBox["\[Tau]", "ii"], "+", 
          SubscriptBox["\[Tau]", "Nee"]}], ")"}]}]}], ")"}]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aie"], " ", 
     SubscriptBox["J", "ei"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "Nee"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "ii"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "Nee"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "ii"], "+", 
          SubscriptBox["\[Tau]", "Nee"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["J", "Nee"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "ei"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "ei"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "Aie"], "+", 
          SubscriptBox["\[Tau]", "ei"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["J", "ii"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "Nee"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "ei"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "Nee"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        RowBox[{"(", 
         RowBox[{
          SubscriptBox["\[Tau]", "ei"], "+", 
          SubscriptBox["\[Tau]", "Nee"], "+", 
          SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "4", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["J", "Nie"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        SubscriptBox["\[Tau]", "ii"], " ", 
        SubscriptBox["\[Tau]", "Nee"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "ii"], " ", 
           SubscriptBox["\[Tau]", "Nee"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "Aie"], " ", 
           RowBox[{"(", 
            RowBox[{
             SubscriptBox["\[Tau]", "ii"], "+", 
             SubscriptBox["\[Tau]", "Nee"]}], ")"}]}]}], ")"}]}]}], ")"}]}], 
    "+", 
    RowBox[{
     SubscriptBox["J", "Aie"], " ", 
     SubscriptBox["J", "ei"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "ii"], " ", 
        SubscriptBox["\[Tau]", "Nee"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "Nee"], " ", 
           SubscriptBox["\[Tau]", "Nie"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "ii"], " ", 
           RowBox[{"(", 
            RowBox[{
             SubscriptBox["\[Tau]", "Nee"], "+", 
             SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], ")"}]}]}],
    ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["J", "Nee"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        SubscriptBox["\[Tau]", "ei"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "ei"], " ", 
           SubscriptBox["\[Tau]", "Nie"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "Aie"], " ", 
           RowBox[{"(", 
            RowBox[{
             SubscriptBox["\[Tau]", "ei"], "+", 
             SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], ")"}]}], 
    "+", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["J", "ii"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "ei"], " ", 
        SubscriptBox["\[Tau]", "Nee"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "Nee"], " ", 
           SubscriptBox["\[Tau]", "Nie"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "ei"], " ", 
           RowBox[{"(", 
            RowBox[{
             SubscriptBox["\[Tau]", "Nee"], "+", 
             SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], 
      ")"}]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "5", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["J", "Nie"], " ", 
     SubscriptBox["\[Tau]", "Aee"], " ", 
     SubscriptBox["\[Tau]", "Aie"], " ", 
     SubscriptBox["\[Tau]", "ii"], " ", 
     SubscriptBox["\[Tau]", "Nee"]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aie"], " ", 
     SubscriptBox["J", "ei"], " ", 
     SubscriptBox["\[Tau]", "Aee"], " ", 
     SubscriptBox["\[Tau]", "ii"], " ", 
     SubscriptBox["\[Tau]", "Nee"], " ", 
     SubscriptBox["\[Tau]", "Nie"]}]}], ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["J", "Nee"], " ", 
     SubscriptBox["\[Tau]", "Aee"], " ", 
     SubscriptBox["\[Tau]", "Aie"], " ", 
     SubscriptBox["\[Tau]", "ei"], " ", 
     SubscriptBox["\[Tau]", "Nie"]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["\[Tau]", "Aie"], " ", 
     SubscriptBox["\[Tau]", "ei"], " ", 
     SubscriptBox["\[Tau]", "Nee"], " ", 
     SubscriptBox["\[Tau]", "Nie"]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "6", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    SubscriptBox["J", "ii"], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       SubscriptBox["\[Tau]", "Aie"], " ", 
       SubscriptBox["\[Tau]", "ei"], " ", 
       SubscriptBox["\[Tau]", "en"], " ", 
       SubscriptBox["\[Tau]", "Nee"], " ", 
       SubscriptBox["\[Tau]", "Nie"]}], "+", 
      RowBox[{
       SubscriptBox["\[Tau]", "Aee"], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          SubscriptBox["\[Tau]", "ei"], " ", 
          SubscriptBox["\[Tau]", "en"], " ", 
          SubscriptBox["\[Tau]", "Nee"], " ", 
          SubscriptBox["\[Tau]", "Nie"]}], "+", 
         RowBox[{
          SubscriptBox["\[Tau]", "Aie"], " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{
             SubscriptBox["\[Tau]", "en"], " ", 
             SubscriptBox["\[Tau]", "Nee"], " ", 
             SubscriptBox["\[Tau]", "Nie"]}], "+", 
            RowBox[{
             SubscriptBox["\[Tau]", "ei"], " ", 
             RowBox[{"(", 
              RowBox[{
               RowBox[{
                SubscriptBox["\[Tau]", "Nee"], " ", 
                SubscriptBox["\[Tau]", "Nie"]}], "+", 
               RowBox[{
                SubscriptBox["\[Tau]", "en"], " ", 
                RowBox[{"(", 
                 RowBox[{
                  SubscriptBox["\[Tau]", "Nee"], "+", 
                  SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], 
           ")"}]}]}], ")"}]}]}], ")"}]}], ">", 
   RowBox[{
    RowBox[{
     SubscriptBox["J", "Nee"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        SubscriptBox["\[Tau]", "ei"], " ", 
        SubscriptBox["\[Tau]", "ii"], " ", 
        SubscriptBox["\[Tau]", "in"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aee"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "ei"], " ", 
           SubscriptBox["\[Tau]", "ii"], " ", 
           SubscriptBox["\[Tau]", "in"], " ", 
           SubscriptBox["\[Tau]", "Nie"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "Aie"], " ", 
           RowBox[{"(", 
            RowBox[{
             RowBox[{
              SubscriptBox["\[Tau]", "ii"], " ", 
              SubscriptBox["\[Tau]", "in"], " ", 
              SubscriptBox["\[Tau]", "Nie"]}], "+", 
             RowBox[{
              SubscriptBox["\[Tau]", "ei"], " ", 
              RowBox[{"(", 
               RowBox[{
                RowBox[{
                 SubscriptBox["\[Tau]", "in"], " ", 
                 SubscriptBox["\[Tau]", "Nie"]}], "+", 
                RowBox[{
                 SubscriptBox["\[Tau]", "ii"], " ", 
                 RowBox[{"(", 
                  RowBox[{
                   SubscriptBox["\[Tau]", "in"], "+", 
                   SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], 
            ")"}]}]}], ")"}]}]}], ")"}]}], "+", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "ei"], " ", 
        SubscriptBox["\[Tau]", "ii"], " ", 
        SubscriptBox["\[Tau]", "in"], " ", 
        SubscriptBox["\[Tau]", "Nee"], " ", 
        SubscriptBox["\[Tau]", "Nie"]}], "+", 
       RowBox[{
        SubscriptBox["\[Tau]", "Aie"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["\[Tau]", "ii"], " ", 
           SubscriptBox["\[Tau]", "in"], " ", 
           SubscriptBox["\[Tau]", "Nee"], " ", 
           SubscriptBox["\[Tau]", "Nie"]}], "+", 
          RowBox[{
           SubscriptBox["\[Tau]", "ei"], " ", 
           RowBox[{"(", 
            RowBox[{
             RowBox[{
              SubscriptBox["\[Tau]", "in"], " ", 
              SubscriptBox["\[Tau]", "Nee"], " ", 
              SubscriptBox["\[Tau]", "Nie"]}], "+", 
             RowBox[{
              SubscriptBox["\[Tau]", "ii"], " ", 
              RowBox[{"(", 
               RowBox[{
                RowBox[{
                 SubscriptBox["\[Tau]", "Nee"], " ", 
                 SubscriptBox["\[Tau]", "Nie"]}], "+", 
                RowBox[{
                 SubscriptBox["\[Tau]", "in"], " ", 
                 RowBox[{"(", 
                  RowBox[{
                   SubscriptBox["\[Tau]", "Nee"], "+", 
                   SubscriptBox["\[Tau]", "Nie"]}], ")"}]}]}], ")"}]}]}], 
            ")"}]}]}], ")"}]}]}], ")"}]}]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "7", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["J", "ii"], " ", 
     SubscriptBox["\[Tau]", "Aee"], " ", 
     SubscriptBox["\[Tau]", "en"], " ", 
     SubscriptBox["\[Tau]", "Nee"]}], "-", 
    RowBox[{
     SubscriptBox["J", "Aee"], " ", 
     SubscriptBox["\[Tau]", "ii"], " ", 
     SubscriptBox["\[Tau]", "in"], " ", 
     SubscriptBox["\[Tau]", "Nee"]}]}], ">", 
   RowBox[{
    SubscriptBox["J", "Nee"], " ", 
    SubscriptBox["\[Tau]", "Aee"], " ", 
    SubscriptBox["\[Tau]", "ii"], " ", 
    SubscriptBox["\[Tau]", "in"]}]}],
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "8", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"],

Cell[BoxData[
 TagBox["True",
  Function[BoxForm`e$, 
   TableForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[BoxData[
 StyleBox[
  RowBox[{
   RowBox[{"Collect", "[", 
    RowBox[{
     RowBox[{"bob2", "[", 
      RowBox[{"[", "9", "]"}], "]"}], ",", "var"}], "]"}], "//", "TableForm"}],
  FontSize->24]], "Input"]
},
ScreenStyleEnvironment->"Presentation",
WindowSize->{1771, 968},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->False,
Deployed->True,
CellContext->Notebook,
TrackCellChangeTimes->False,
FrontEndVersion->"9.0 for Linux x86 (64-bit) (February 7, 2013)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1485, 35, 210, 7, 52, "Input"],
Cell[1698, 44, 472, 18, 78, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2207, 67, 210, 7, 52, "Input"],
Cell[2420, 76, 1344, 42, 116, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3801, 123, 210, 7, 52, "Input"],
Cell[4014, 132, 2902, 90, 116, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6953, 227, 210, 7, 52, "Input"],
Cell[7166, 236, 3079, 98, 116, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10282, 339, 210, 7, 52, "Input"],
Cell[10495, 348, 1140, 34, 78, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[11672, 387, 210, 7, 52, "Input"],
Cell[11885, 396, 4560, 128, 153, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16482, 529, 210, 7, 52, "Input"],
Cell[16695, 538, 629, 20, 78, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17361, 563, 210, 7, 52, "Input"],
Cell[17574, 572, 92, 3, 78, "Output"]
}, Open  ]],
Cell[17681, 578, 210, 7, 52, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature XvDVV@Z0nAGIxAKuWcjswoIR *)
