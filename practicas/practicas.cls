% License missing **********************************************
%
% Repository:  **************************************
%
% Notes
% ~~~~~
% The following tags are used:
% HW        :   User-level class commands
% CeHW@     :   The generic tag used through the class
% CeHWl@    :   The tag used to define labels and translations (except for the list of problems that follows the example of figure)
% CeHWf@    :   The tag used to define the fonts styles
%

\ProvidesClass{practicas}[2022/01/15 V3.0 by Amado Cabrera and Mariana Pérez]
\NeedsTeXFormat{LaTeX2e}

\typeout{-- .}
\typeout{-- .}%***************************************

%% Class version %%
% CabesHW.cls version numbers, provided as of V1.3
% These values serve as a way a .tex file can
% determine if the new features are provided.
% The version number of this CabesHW.cls can be obtained from 
% these values. i.e., V1.4
\def\HWvMajor{3}
\def\HWvMinor{1}

%% Class option conditionals %%
\newif\ifCeHW@table     \CeHW@tablefalse   % Table packages
\newif\ifCeHW@theorem   \CeHW@theoremfalse % To load my theorems
\newif\ifCeHW@classicMT \CeHW@classicMTfalse % Classic Maketitle

%% Declaring and loading the class options %%
\DeclareOption{Ce-table}{\CeHW@tabletrue}
\DeclareOption{Ce-theorem}{\CeHW@theoremtrue}
\DeclareOption{Ce-classicMT}{\CeHW@classicMTtrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*\relax
\LoadClass{article}

%% Loading packages %%
\RequirePackage{iftex} % Manage the loading depending on the engine

\ifPDFTeX % PDFTeX -> (latex and pdflatex)
    \RequirePackage[T1]{fontenc} % Base for good pdf encoding
    \RequirePackage{lmodern}
\fi\ifLuaTeX % LuaTeX -> (and extensions such as LuaHBTeX)
    \RequirePackage{fontspec} % For different font usage
\fi
\RequirePackage{geometry} % Better layout
\RequirePackage{xcolor} % Better design
\RequirePackage{microtype} % Improve micro-typography
\RequirePackage{graphicx} % Images with keyval support
\RequirePackage{float} % Non floating figures, etc
\RequirePackage[shortlabels,inline]{enumitem} % Better lists
\RequirePackage{titlesec} % Better control on sectioning commands
\RequirePackage{titletoc} % Better control on the TOC
\RequirePackage{amsmath} % Math
\RequirePackage{amsthm} % More Math
\RequirePackage{amssymb} % Math fonts
\RequirePackage{dsfont} % Load the Dsfont for numeric sets
\RequirePackage{mathrsfs} % Load the rsfs

% Conditional loadings
\if@twocolumn
    \relax
\else
    \RequirePackage{multicol} % Easy multiple columns
\fi
\ifCeHW@table
    \RequirePackage{array} % Better control on table columns
\fi


%% Packages setup %%
% geometry package setup
\geometry{
    letterpaper,
    % showframe,
    margin=1in
}
% babel package setup
\newcommand\CeHW@iflangp[1]{%
  \@ifpackageloaded{babel}
    {#1}
    {%
      \@ifpackageloaded{polyglossia}
        {#1}
        {}%
    }%
}

% Default
\def\CeHWl@author{Name:\space}
\def\CeHWl@date{}
\def\CeHWl@institute{}
\def\CeHWl@course{Course:\space}
\def\CeHWl@id{Id:\space}
\def\CeHWl@professor{Professor:\space}
\def\CeHWl@email{Email:\space}

\AtBeginDocument{\CeHW@iflangp{
% English
\addto\captionsenglish{
\gdef\CeHWl@author{Name:\space}
\gdef\CeHWl@date{}
\gdef\CeHWl@institute{}
\gdef\CeHWl@course{}
\gdef\CeHWl@id{Id:\space}
\gdef\CeHWl@professor{Supervisado por:\space}
\gdef\CeHWl@email{E--mail:\space}
\gdef\listCeHWproblemname{List of Problems}
\gdef\CeHWproblemname{Problem}
}
% Spanish
\addto\captionsspanish{
\gdef\CeHWl@author{Nombre:\space}
\gdef\CeHWl@date{}
\gdef\CeHWl@institute{}
\gdef\CeHWl@course{}
\gdef\CeHWl@id{Carné:\space}
\gdef\CeHWl@professor{Asesorado por:\space}
\gdef\CeHWl@email{E--mail:\space}
\gdef\listCeHWproblemname{Índice de problemas}
\gdef\CeHWproblemname{Problema}
}
}}

% xcolor package setup
\let\CeHW@oldpagecolor\pagecolor

\gdef\CeHW@thepagecolor{white}
\renewcommand{\pagecolor}[1]{\gdef\CeHW@thepagecolor{#1}%
\CeHW@oldpagecolor{#1}}

\definecolor{HWcolor-1}{HTML}{000000}%
\definecolor{HWcolor-2}{HTML}{000000}%
\definecolor{HWcolor-3}{HTML}{000000}%
\definecolor{HWcolor-4}{HTML}{000000}%
\definecolor{HWcolor-5}{HTML}{000000}%
\def\CeHW@setPalette#1,#2,#3,#4,#5\ENDCMD{%
    \definecolor{HWcolor-1}{HTML}{#1}%
    \definecolor{HWcolor-2}{HTML}{#2}%
    \definecolor{HWcolor-3}{HTML}{#3}%
    \definecolor{HWcolor-4}{HTML}{#4}%
    \definecolor{HWcolor-5}{HTML}{#5}%
}
\newcommand{\setPalette}[1]{\CeHW@setPalette#1\ENDCMD}


% enumitem package setup
% star items
% Code from egreg (https://tex.stackexchange.com/a/167517/234538) with some modifications to match the programming style of the class
\newcommand{\staritem}{\global\CeHW@startrue\item}
\newcommand{\CeHW@perhapsstar}{%
  \ifCeHW@star\raisebox{0.1em}{$\star$}\global\CeHW@starfalse\fi
}
\newif\ifCeHW@star
% enumerate config
\let\CeHW@enumi\labelenumi
\let\CeHW@enumii\labelenumii
\let\CeHW@enumiii\labelenumiii
\let\CeHW@enumiv\labelenumiv
\setlist[enumerate,1]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@enumi}%
}
\setlist[enumerate,2]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@enumii}%
}
\setlist[enumerate,3]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@enumiii}%
}
\setlist[enumerate,4]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@enumiv}%
}
% itemize config
\let\CeHW@itemi\labelitemi
\let\CeHW@itemii\labelitemii
\let\CeHW@itemiii\labelitemiii
\let\CeHW@itemiv\labelitemiv
\setlist[itemize,1]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@itemi}%
}
\setlist[itemize,2]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@itemii}%
}
\setlist[itemize,3]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@itemiii}%
}
\setlist[itemize,4]{%
    label={\protect\CeHW@perhapsstar\,\CeHW@itemiv}%
}
% description config
\setlist[description,1]{font={}}
\setlist[description,2]{font={}}
\setlist[description,3]{font={}}
\setlist[description,4]{font={}}


% titlesec package setup
\def\@part[#1]#2{%                  % exception with part
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont\color{HWcolor-1}
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname\nobreakspace\thepart
       \par\nobreak
     \fi
     \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont\color{HWcolor-1}
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}

\titleformat{\section}%
    {\normalfont\Large\bfseries\color{HWcolor-2}}%
    {\thesection}{1em}{}
\titlespacing*{\section}{0pt}%
    {3ex plus 1ex minus .2ex}%
    {2.3ex plus .2ex}

\titleformat{\subsection}%
    {\normalfont\large\bfseries\color{HWcolor-3}}%
    {\thesubsection}{1em}{}
\titlespacing*{\subsection}{0pt}%
    {2.75ex plus 1ex minus .2ex}%
    {1.5ex plus .2ex}

\titleformat{\subsubsection}%
    {\normalfont\normalsize\bfseries\color{HWcolor-4}}%
    {\thesubsubsection}{1em}{}
\titlespacing*{\subsubsection}{0pt}%
    {2.75ex plus 1ex minus .2ex}%
    {1.5ex plus .2ex}

\titleformat{\paragraph}[runin]%
    {\normalfont\normalsize\bfseries\color{HWcolor-5}}%
    {\theparagraph}{1em}{}
\titlespacing*{\paragraph}{0pt}%
    {3.25ex plus 1ex minus .2ex}%
    {1em}

\let\subparagraph\relax

% titlesec package setup
\titlecontents{part}[0pt]
    {\large\bfseries\bigskip}
    {}
    {}
    {\hfill\contentspage}

\titlecontents{section}[0pt]
    {\bfseries\medskip}
    {\hspace{1.3em}\contentslabel{1.3em}}
    {}
    {\hfill\contentspage}

\titlecontents{subsection}[1.3em]
    {\smallskip}
    {\hspace{2.3em}\contentslabel{2.3em}}
    {}
    {\titlerule*[0.5pc]{.}\contentspage}

\titlecontents{subsubsection}[3.6em]
    {}
    {\hspace{3.2em}\contentslabel{3.2em}}
    {}
    {\titlerule*[0.5pc]{.}\contentspage}

\titlecontents{figure}[1.5em]
    {}
    {\hspace{2.3em}\contentslabel{2.3em}}
    {}
    {\titlerule*[0.5pc]{.}\contentspage}

\titlecontents{table}[1.5em]
    {}
    {\hspace{2.3em}\contentslabel{2.3em}}
    {}
    {\titlerule*[0.5pc]{.}\contentspage}

\let\TOC\tableofcontents

\renewcommand{\tableofcontents}{
\titlespacing*{\section}{0pt}%
    {3ex plus 1ex minus .2ex}%
    {0.5ex plus .1ex}%
\TOC%
\titlespacing*{\section}{0pt}%
    {3ex plus 1ex minus .2ex}%
    {2.3ex plus .2ex}}


% amsmath and mathematical commands
% Extra math operators
    % trigonometric
\DeclareMathOperator{\arccosh}{arccosh}
\DeclareMathOperator{\arcsinh}{arcsinh}
    % matrices
\DeclareMathOperator{\Tr}{Tr}
    % vectors
\let\division\div
\let\div\relax
\DeclareMathOperator{\grad}{\nabla}
\DeclareMathOperator{\laplacian}{\nabla^2}
\DeclareMathOperator{\div}{\nabla\cdot}
\DeclareMathOperator{\curl}{\nabla\times}

% named equations
\newlength{\CeHW@nequation}
\setlength{\CeHW@nequation}{0.25\columnwidth}
\def\nequationfont{\bfseries}
\newenvironment{nequation}[1]{%
\begin{list}{}{%
    \topsep=\abovedisplayskip\relax%
    \partopsep=\z@\relax%
    \parsep=\z@\relax%
    \rightmargin=\z@\relax%
    \leftmargin=\z@\relax%
    \listparindent=\z@\relax%
    \itemsep=\z@\relax%
}
\item\begin{minipage}{\CeHW@nequation}\nequationfont%
#1%
\end{minipage}%
\hspace{-\CeHW@nequation}%
\begin{minipage}{\columnwidth}%
\begin{equation}
}{%
\end{equation}
\end{minipage}%
\end{list}
}
\newenvironment{nequation*}[1]{%
\begin{list}{}{%
    \topsep=\abovedisplayskip\relax%
    \partopsep=\z@\relax%
    \parsep=\z@\relax%
    \rightmargin=\z@\relax%
    \leftmargin=\z@\relax%
    \listparindent=\z@\relax%
    \itemsep=\z@\relax%
}
\item\begin{minipage}{\CeHW@nequation}\nequationfont%
#1%
\end{minipage}%
\hspace{-\CeHW@nequation}%
\begin{minipage}{\columnwidth}%
\begin{equation*}
}{%
\end{equation*}
\end{minipage}%
\end{list}
}

% diferentials and derivatives
\newcommand{\dif}{\,\mathrm{d}}
\newcommand{\pdif}{\,\partial}

\def\td{\@ifstar\CeHW@tdI\CeHW@tdD}
\newcommand{\CeHW@tdD}[3][]{%
    \frac{\mathrm{d}^{#1}#2}{\mathrm{d}#3^{#1}}%
}
\newcommand{\CeHW@tdI}[3][]{%
    \mathrm{d}^{#1}#2/\mathrm{d}#3^{#1}%
}

\def\pd{\@ifstar\CeHW@pdI\CeHW@pdD}
\newcommand{\CeHW@pdD}[3][]{%
    \frac{\partial^{#1}#2}{\partial\CeHW@pdDenominator#3,\CeHW@delim}%
}
\newcommand{\CeHW@pdI}[3][]{%
    \partial^{#1}#2/\partial\CeHW@pdDenominator#3,\CeHW@delim %
}
\def\CeHW@pdDenominator#1,#2\CeHW@delim{%
    \ifx\CeHW@undefined#2\CeHW@undefined%
        {#1}%
    \else%
        {#1\,\partial\CeHW@pdDenominator#2\CeHW@delim}%
    \fi%
}

% Bra-Ket notation
\DeclareRobustCommand{\bra}{\@ifstar\CeHW@braI\CeHW@braD}
\newcommand{\CeHW@braI}[1]{%
    \@ifnextchar\ket{\CeHW@braIi{#1}}{\left\langle\smash{#1}\right\rvert}%
}
\def\CeHW@braIi#1\ket{\@ifnextchar*{\CeHW@braII{#1}}{\CeHW@braketIi{#1}}}
\def\CeHW@braII#1*#2{\CeHW@braketIi{#1}{#2}}
\newcommand{\CeHW@braD}[1]{%
    \@ifnextchar\ket{\CeHW@braDi{#1}}{\left\langle{#1}\right\rvert}%
}
\def\CeHW@braDi#1\ket{\@ifnextchar*{\CeHW@braDI{#1}}{\CeHW@braketDi{#1}}}
\def\CeHW@braDI#1*#2{\CeHW@braketIi{#1}{#2}}

\DeclareRobustCommand{\ket}{\@ifstar\CeHW@ketI\CeHW@ketD}
\newcommand{\CeHW@ketI}[1]{%
    \left\lvert\smash{#1}\right\rangle%
}
\newcommand{\CeHW@ketD}[1]{%
    \left\lvert{#1}\right\rangle%
}

\newcommand{\braket}{\@ifstar\CeHW@braketI\CeHW@braketD}
\newcommand{\CeHW@braketI}[1]{%
    \@ifnextchar\bgroup{\CeHW@braketIi{#1}}{\CeHW@braketIi{#1}{#1}}%
}
\newcommand{\CeHW@braketIi}[2]{%
    \left\langle\smash{#1}\middle\vert\smash{#2}\right\rangle%
}
\newcommand{\CeHW@braketD}[1]{%
    \@ifnextchar\bgroup{\CeHW@braketDi{#1}}{\CeHW@braketDi{#1}{#1}}%
}
\newcommand{\CeHW@braketDi}[2]{%
    \left\langle{#1}\middle\vert{#2}\right\rangle%
}

\newcommand{\dyad}{\@ifstar\CeHW@dyadI\CeHW@dyadD}
\newcommand{\CeHW@dyadI}[1]{%
    \@ifnextchar\bgroup{\CeHW@dyadIi{#1}}{\CeHW@dyadIi{#1}{#1}}%
}
\newcommand{\CeHW@dyadIi}[2]{%
    \left\lvert\smash{#1}\middle\rangle\!\middle\langle\smash{#2}\right\rvert%
}
\newcommand{\CeHW@dyadD}[1]{%
    \@ifnextchar\bgroup{\CeHW@dyadDi{#1}}{\CeHW@dyadDi{#1}{#1}}%
}
\newcommand{\CeHW@dyadDi}[2]{%
    \left\lvert{#1}\middle\rangle\!\middle\langle{#2}\right\rvert%
}
\let\ketbra\dyad

\newcommand{\expval}{\@ifstar\CeHW@expvalI\CeHW@expvalD}
\newcommand{\CeHW@expvalI}[1]{%
    \@ifnextchar\bgroup{\CeHW@expvalIi{#1}}{\left\langle\smash{#1}\right\rangle}%
}
\newcommand{\CeHW@expvalIi}[2]{%
    \left\langle\smash{#2}\middle\vert\smash{#1}\middle\vert\smash{#2}\right\rangle%
}
\newcommand{\CeHW@expvalD}[1]{%
    \@ifnextchar\bgroup{\CeHW@expvalDi{#1}}{\left\langle{#1}\right\rangle}%
}
\newcommand{\CeHW@expvalDi}[2]{%
    \left\langle{#2}\middle\vert{#1}\middle\vert{#2}\right\rangle%
}

\newcommand{\mel}{\@ifstar\CeHW@melI\CeHW@melD}
\newcommand{\CeHW@melI}[3]{%
    \left\langle\smash{#1}\middle\vert\smash{#2}\middle\vert\smash{#3}\right\rangle%
}
\newcommand{\CeHW@melD}[3]{%
    \left\langle{#1}\middle\vert{#2}\middle\vert{#3}\right\rangle%
}

% Matrices
\newcommand{\mat}[2][p]{%
    \begin{#1matrix}%
    #2%
    \end{#1matrix}%
}

\newcommand{\imat}[2][p]{%
\count\@ne=0
\toks@{}
\loop
    \count@=0
    \advance\count\@ne by 1\relax
    {%
    \loop
        \advance\count@ by 1\relax
        \edef\addtoks{%
            \ifnum\count@=\count\@ne%
                1%
            \else%
                0%
            \fi%
            \ifnum\count@=#2%
            \else%
                &%
            \fi%
        }
        \global\toks@ = \expandafter{\the\expandafter\toks@\addtoks}
    \ifnum\count@<#2\repeat
    }
    \edef\addtoks{%
        \ifnum\count@=#2%
            \else%
                \cr%
            \fi%
    }
    \toks@ = \expandafter{\the\expandafter\toks@\addtoks}
\ifnum\count\@ne<#2\repeat
\begin{#1matrix}%
\the\toks@%
\end{#1matrix}%
}

\newcommand{\xmat}[4][p]{%
\count\@ne=0
\toks@{}
\loop
    \count@=0
    \advance\count\@ne by 1\relax
    {%
    \loop
        \advance\count@ by 1\relax
        \edef\addtoks{%
            #2
            \ifnum\count@=#4%
            \else%
                &%
            \fi%
        }
        \global\toks@ = \expandafter{\the\expandafter\toks@\addtoks}
    \ifnum\count@<#4\repeat
    }
    \edef\addtoks{%
        \ifnum\count@=#3%
            \else%
                \cr%
            \fi%
    }
    \toks@ = \expandafter{\the\expandafter\toks@\addtoks}
\ifnum\count\@ne<#3\repeat
\begin{#1matrix}%
\the\toks@%
\end{#1matrix}%
}

\newcommand{\zmat}[3][p]{\xmat[#1]{0}{#2}{#3}}

\newcounter{nummat}%
\newcommand{\paulimat}[2][p]{%
\setcounter{nummat}{#2}%
\ifcase\thenummat%
    \edef\@val{1&0\\0&1}
\or
    \edef\@val{0&1\\1&0}
\or
    \edef\@val{0&-i\\i&0}
\or
    \edef\@val{1&0\\0&-1}
\fi
\begin{#1matrix}%
\@val
\end{#1matrix}%
}

\newcommand{\dmat}[2][p]{%
\count\@ne=0
\toks@{}
\@for\i@val:=#2\do{%
    \count@=0
    {
    \loop
        \edef\addtoks{\ifnum\count\@ne=0\else&\fi}
        \global\toks@ = \expandafter{\the\expandafter\toks@\addtoks}
        \advance\count@ by 1\relax%
    \ifnum\count@<\count\@ne\repeat
    }
    \toks@ = \expandafter{\the\expandafter\toks@\i@val\cr}
    \advance\count\@ne by 1\relax
}
\begin{#1matrix}
\the\toks@
\end{#1matrix}%
}

% Math boxes
\newcommand{\Cboxed}[2][HWcolor-3]{%
    \fcolorbox{#1}{\CeHW@thepagecolor}{\m@th$\displaystyle#2$}%
}
% Code from campa (https://tex.stackexchange.com/a/610299/234538) with some modifications to match the programming style of the class
\newcommand*\CAboxed[2][HWcolor-3]{%
   \let\bgroup{\romannumeral-`}%
   \@CAboxed{#1}#2&&\ENDDNE
}

\def\@CAboxed#1#2&#3&#4\ENDDNE{%
  \ifnum0=`{}\fi \setbox \z@
    \hbox{$\displaystyle#2{}\m@th$\kern\fboxsep \kern\fboxrule}%
    \edef\@tempa {\kern \wd\z@ &\kern -\the\wd\z@ \fboxsep 
        \the\fboxsep \fboxrule\the\fboxrule }\@tempa \Cboxed[#1] {#2#3}%
}


% amsthm configuration
% styles
\newtheoremstyle{HWbreak}
    {}
    {}
    {\slshape}
    {\z@}
    {\bfseries}
    {}
    {\z@}
    {\thmname{#1}\thmnumber{~#2}:\thmnote{ \mdseries#3.}\\[2\p@]}
\newtheoremstyle{HWbreak*}
    {}
    {}
    {\slshape}
    {\z@}
    {\bfseries}
    {}
    {\z@}
    {\ifx\HW@undefined#3\HW@undefined%
        {\thmname{#1}}%
    \else%
        {\thmname{#3}}%
    \fi:\\[2\p@]}
\renewcommand{\qedsymbol}{$\mathcal{Q.E.D.}$}

% environments
\def\CeHWl@theorem{Teorema}
\def\CeHWl@definition{Definición}
\def\CeHWl@corollary{Corolario}
\def\CeHWl@lemma{Lemma}

\ifCeHW@theorem
\theoremstyle{HWbreak}
\newtheorem{theorem}{\CeHWl@theorem}
\newtheorem{definition}{\CeHWl@definition}
\newtheorem{corollary}{\CeHWl@corollary}
\newtheorem{lemma}{\CeHWl@lemma}
\theoremstyle{HWbreak*}
\newtheorem*{theorem*}{\CeHWl@theorem}
\theoremstyle{plain}
\fi

%% Problems stuff %%
\newcounter{CeHWproblem}
\renewcommand\theCeHWproblem{\@arabic\c@CeHWproblem}
\def\ext@CeHWproblem{Celop}
\def\fnum@CeHWproblem{\CeHWproblemname\nobreakspace\theCeHWproblem}
\newcommand\listCeHWproblemname{List of Problems}
\newcommand\CeHWproblemname{Problem}
\newcommand{\CeHW@problemformat}[1]{\nobreakspace(#1)}
\def\ttll@CeHWproblem{1}

\newcommand\listofproblems{%
    \section*{\listCeHWproblemname}%
    \@mkboth{\MakeUppercase\listCeHWproblemname}%
        {\MakeUppercase\listCeHWproblemname}%
    \@starttoc{Celop}%
}

\titlecontents{CeHWproblem}[1.5em]
    {}
    {\hspace{2.3em}\contentslabel{2.3em}}
    {}
    {\titlerule*[0.5pc]{.}\contentspage}

\newenvironment{prob}[1][]{%
    \begin{trivlist}\item
    \titlespacing*{\section}{0pt}%
        {0pt}%
        {2.3ex plus .2ex}
    \vskip 3ex plus 1ex minus .2ex\relax
    \refstepcounter{CeHWproblem}%
    \addcontentsline{\ext@CeHWproblem}{CeHWproblem}%
        {\protect\numberline{\theCeHWproblem}{\ignorespaces\fnum@CeHWproblem%
            \ifx\CeHW@undefined#1\CeHW@undefined%
            \else%
                \CeHW@problemformat{#1}%
            \fi%
        }}%
    \section*{%
        \fnum@CeHWproblem%
        \ifx\CeHW@undefined#1\CeHW@undefined%
        \else%
            \CeHW@problemformat{#1}%
        \fi%
    }%
}{%
    \end{trivlist}%
}%

\newenvironment{prob*}[1][]{%
    \begin{trivlist}\item%
    \ifx\CeHW@undefined#1\CeHW@undefined%
    \else%
        \section*{#1}%
    \fi%
}{%
    \end{trivlist}%
}%

\newlist{ListProbListType}{enumerate}{1}
\setlist[ListProbListType,1]{label={\alph*)}}

\newenvironment{listprob}[1][]{%
    \renewenvironment{ans}[1][0pt]{%
    \begin{HWadjustwidth}{##1}{##1}}%
    {\end{HWadjustwidth}}
    %
    \begin{trivlist}\item
    \titlespacing*{\section}{0pt}%
        {0pt}%
        {2.3ex plus .2ex}
    \vskip 3ex plus 1ex minus .2ex\relax
    \refstepcounter{CeHWproblem}%
    \addcontentsline{\ext@CeHWproblem}{CeHWproblem}%
        {\protect\numberline{\theCeHWproblem}{\ignorespaces\fnum@CeHWproblem%
            \ifx\CeHW@undefined#1\CeHW@undefined%
            \else%
                \CeHW@problemformat{#1}%
            \fi%
        }}%
    \section*{%
        \fnum@CeHWproblem%
        \ifx\CeHW@undefined#1\CeHW@undefined%
        \else%
            \CeHW@problemformat{#1}%
        \fi%
    }%
    \begin{ListProbListType}
}{%
    \end{ListProbListType}
    \end{trivlist}%
}%

\newenvironment{listprob*}[1][]{%
    \renewenvironment{ans}[1][0pt]{%
    \begin{HWadjustwidth}{##1}{##1}}%
    {\end{HWadjustwidth}}
    %
    \begin{trivlist}\item%
    \ifx\CeHW@undefined#1\CeHW@undefined%
    \else%
        \section*{#1}%
    \fi%
    \begin{ListProbListType}
}{%
    \end{ListProbListType}
    \end{trivlist}%
}%

\newenvironment{HWadjustwidth}[2]{%
    \begin{list}{}{%
    \topsep\z@%
    \listparindent\parindent%
    \parsep\parskip%
    \ifx\HW@undefined#1\HW@undefined%
        \setlength{\leftmargin}{\z@}%
    \else%
        \setlength{\leftmargin}{#1}%
    \fi%
    \ifx\HW@undefined#2\HW@undefined%
        \setlength{\rightmargin}{\z@}%
    \else%
        \setlength{\rightmargin}{#2}%
    \fi%
    }
    \item[]}%
    {\end{list}}

\newenvironment{ans}[1][1.5em]{%
    \begin{HWadjustwidth}{#1}{#1}}%
    {\end{HWadjustwidth}}

%% Maketitle stuff %%
% Auxiliary commands
\newcommand{\CeHW@isUndef}[2]{%
\ifx#1\HW@undefined%
    {}%
\else%
    {#2}%
\fi%
}
\newcommand{\CeHW@maketitletype}[2][]{{\csname CeHWl@#2\endcsname}{#1\csname CeHWf@#2\endcsname\csname CeHW@#2\endcsname}}

\def\CeHWf@author{\itshape}
\def\CeHWf@date{}

\newcommand{\institute}[1]{\gdef\CeHW@institute{#1}}
\def\CeHWf@institute{\scshape}

\newcommand{\course}[1]{\gdef\CeHW@course{#1}}
\def\CeHWf@course{\itshape}

\newcommand{\id}[1]{\gdef\CeHW@id{#1}}
\def\CeHWf@id{}

\newcommand{\professor}[1]{\gdef\CeHW@professor{#1}}
\def\CeHWf@professor{\itshape}

\newcommand{\email}[1]{\gdef\CeHW@email{#1}}
\def\CeHWf@email{}

\renewcommand{\and}{\&\space}

% Maketitle commands
\ifCeHW@classicMT
\else
\if@titlepage
\gdef\CeHWl@institute{}
\gdef\CeHWl@course{}
\gdef\CeHWl@id{}
\gdef\CeHWl@professor{}
\gdef\CeHWl@email{}
  \renewcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \begin{center}%
    \CeHW@isUndef{\CeHW@institute}{\Large\CeHW@maketitletype{institute}}%
    \null\vfill%
    \begin{minipage}{2\p@}%
    \vrule width 2\p@ height 0.85\textheight\relax%
    \end{minipage}\@killglue\kern 20\p@\relax%
    \begin{minipage}{0.65\textwidth}%
    {\Huge \@title}%
    \vskip 3em%
    {\Large\CeHWf@author\@author}%
    \vskip 1em%
    \CeHW@isUndef{\CeHW@id}{\Large\CeHW@maketitletype{id}}%
    \vskip 2em%
    {\large \@date}% 
    \vskip 0.425\textheight\relax\null
    \CeHW@isUndef{\CeHW@course}{\large\CeHW@maketitletype{course}}%
    \vskip 1em%
    \CeHW@isUndef{\CeHW@professor}{\large\CeHW@maketitletype{professor}}%
    \end{minipage}\kern 22\p@\relax%
    \vfill\null%
  \end{center}
  \@thanks%
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax%
  \global\let\maketitle\relax%
  \global\let\@thanks\@empty%
  \global\let\@author\@empty%
  \global\let\@date\@empty%
  \global\let\@title\@empty%
  \global\let\title\relax%
  \global\let\author\relax%
  \global\let\date\relax%
  \global\let\and\relax%
}
\else
\renewcommand\maketitle{\par%
  \begingroup%
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage%
  \begin{center}%
  \let \footnote \thanks%
    \CeHW@isUndef{\CeHW@institute}{\CeHW@maketitletype{institute}}
    \vskip 5\p@\relax%
    \hrule width \textwidth\relax%
    \vskip 10\p@\relax%
    {\huge \@title \\[2\p@]}%
    \vskip 1em%
    \CeHWl@author{\CeHWf@author\@author}\\[\p@]%
    \vskip 2.5\p@\relax
    \CeHW@isUndef{\CeHW@id}{\CeHW@maketitletype{id}}%
    \vskip 8\p@\relax%
    \hrule width \textwidth\relax%
    \vskip 7\p@\relax%
    \CeHW@isUndef{\CeHW@course}{\CeHW@maketitletype{course}%
    \CeHW@isUndef{\CeHW@professor}{\space--\space}}%
    \CeHW@isUndef{\CeHW@professor}{\CeHW@maketitletype{professor}}%
    \vskip .3em%
    {\CeHWf@date\@date}%
  \end{center}%
  \vskip 1.5em}
\fi
\fi

\endinput
