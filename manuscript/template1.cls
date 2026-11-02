\NeedsTeXFormat{LaTeX2e} \ProvidesClass{template}
\RequirePackage{setspace}

\LoadClass[11pt, a4paper]{report}

\oddsidemargin 14.6truemm \evensidemargin 0truemm
\marginparwidth 40pt \marginparsep 10pt
\topmargin 0truemm \headsep 10truemm
\textheight 224truemm \textwidth 145truemm

\brokenpenalty=10000

\onehalfspacing



\usepackage{amsmath}
\usepackage{mfirstuc}
\usepackage{graphicx}
%\usepackage{geometry}
\usepackage[bindingoffset=10mm,bottom=4cm]{geometry}
%\usepackage{float}
\usepackage{booktabs}
\usepackage[explicit]{titlesec}
\usepackage[hidelinks]{hyperref}
\usepackage[nameinlink]{cleveref}
\usepackage{pgfplotstable}
\pgfplotsset{compat=1.18}
\usepackage[per-mode=symbol]{siunitx}
\usepackage[toc,title]{appendix}
\usepackage{url}
\usepackage[section]{placeins}
\usepackage{nomencl}
\usepackage{floatrow}
\floatsetup[figure]{capposition=bottom}
\floatsetup[table]{capposition=top}
\floatplacement{figure}{H}
\floatplacement{table}{H}



\makenomenclature

\AtBeginDocument{\maketitle}
\bibliographystyle{unsrt}

\newcommand\abst[1]{\newcommand\zzabst{#1}}

\newcommand*{\project}[1]{\gdef\@project{#1}}
\newcommand*{\@project}{Graduation Project I}
\newcommand*{\supervisor}[1]{\gdef\@supervisor{#1}}
\newcommand*{\myabstract}[1]{\gdef\@myabstract{#1}}
\newcommand*{\faculty}[1]{\gdef\@faculty{#1}}
\newcommand*{\@faculty}{Engineering}
\newcommand*{\department}[1]{\gdef\@department{#1}}
\newcommand*{\@department}{Mechanical Engineering}

\makeatletter
\AtEndEnvironment{figure}{\gdef\there@is@a@figure{}} 
\AtEndDocument{\ifdefined\there@is@a@figure\label{fig:was:used:in:doc}\fi} 
\newcommand{\conditionalLoF}{\@ifundefined{r@fig:was:used:in:doc}{}{\newpage\listoffigures}}%
\makeatother

\makeatletter
\AtEndEnvironment{table}{\gdef\there@is@a@table{}} 
\AtEndDocument{\ifdefined\there@is@a@table\label{tab:was:used:in:doc}\fi} 
\newcommand{\conditionalLoT}{\@ifundefined{r@tab:was:used:in:doc}{}{\newpage\listoftables}}%
\makeatother

\makeatletter
\renewcommand*{\maketitle}{%
  \begin{titlepage}
  \newgeometry{bindingoffset=10mm,bottom=3cm}
  {\centering\includegraphics[scale=0.2]{logo}\par}
  \vspace{1cm}
  \centering
  {\scshape\LARGE S.N. Bose National Centre for Basic Sciences \par}
  \vspace{1cm}
  {\Large\xcapitalisewords Saltlake, Kolkata-700106 \unskip\strut\par}
  \vspace{1cm}
  {\Large\xcapitalisewords Department of {\@department}\unskip\strut\par}
  \vspace{2cm}
  {\scshape\Large\@project\unskip\strut\par}
  \vspace{1.5cm}
  {\huge\bfseries\xcapitalisewords{\@title}\unskip\strut\par}
  \vspace{2cm}
  {\Large\xcapitalisewords{\@author}\unskip\strut\par}  
  \vfill
  Supervised by\par
  Dr. Suman Chakrabarty \unskip\strut\par
  \vfill
  {\large \@date\par}
  \end{titlepage}
  \newpage
  \pagenumbering{roman}
  
\begin{abstract}
\zzabst
\end{abstract}

\tableofcontents
%\listoffigures
%\listoftables
\conditionalLoF
\conditionalLoT
\printnomenclature
\newpage
\pagenumbering{arabic}
}
\makeatother

\usepackage{etoolbox}
\patchcmd{\abstract}{\null\vfil}{}{}{}

\renewcommand{\bibname}{\raggedright References}

\titleformat{\chapter}[display]{\bfseries\centering}{\huge \chaptertitlename\ \thechapter}{1em}{\Huge #1}
\renewcommand{\contentsname}{\raggedright Table of Contents}
\renewcommand{\listfigurename}{\raggedright List of Figures}
\renewcommand{\listtablename}{\raggedright List of Tables}

\pgfplotstableset{
    every head row/.style={
    before row=\toprule,after row=\midrule},
    every last row/.style={
    after row=\bottomrule},
    col sep = &,
    row sep=\\,
    string type,
}

\newcommand{\entries}{\pgfplotstabletypeset}

\newtoggle{appendixbib}
\pretocmd\appendix{\toggletrue{appendixbib}\bibliography{ref}}{}{}
\AtEndDocument{\iftoggle{appendixbib}{}{\bibliography{ref}}}

%\let\origfigure\figure
%\let\endorigfigure\endfigure

%\renewenvironment{figure}[1][H]{%
%	\origfigure[#1]%
%	\centering
%}{%
%	\endorigfigure
%}

%\let\origtable\table
%\let\endorigtable\endtable

%\renewenvironment{table}[1][H]{%
%	\origtable[#1]%
%	\centering
%}{%
%	\endorigtable
%}

