#!/usr/bin/env python3
"""Build the Chinese method guide with XeLaTeX and real math environments.

The Markdown source remains the readable content source. ``[eq]`` lines are
converted through an explicit map rather than emitted as plain text.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


FORMULAS = {
    "P_eq : Ω_o(P,T) = Ω_c(P,T)": r"P_{\mathrm{eq}}:\quad \Omega_o(P,T)=\Omega_c(P,T)",
    "E_i[kJ/mol-f.u.] = 96.4853321233 × E_i[eV/CIF cell] / Z_i": r"E_i[\mathrm{kJ\,mol^{-1}\,f.u.^{-1}}]=96.4853321233\,\frac{E_i[\mathrm{eV\,CIF\ cell}]}{Z_i}",
    "U_ads = U_host–guest + U_guest–guest": r"U_{\mathrm{ads}}=U_{\mathrm{host-guest}}+U_{\mathrm{guest-guest}}",
    "U_LJ(r) = 4ε[(σ/r)^12 − (σ/r)^6]": r"U_{\mathrm{LJ}}(r)=4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]",
    "ΔE_static = E_o/Z_o − E_c/Z_c": r"\Delta E_{\mathrm{static}}=\frac{E_o}{Z_o}-\frac{E_c}{Z_c}",
    "K_ij = k_parallel [ n n^T + α(I − n n^T) ]": r"\mathbf K_{ij}=k_{\parallel}\left[\mathbf n\mathbf n^{\mathsf T}+\alpha\left(\mathbf I-\mathbf n\mathbf n^{\mathsf T}\right)\right]",
    "z_vib = exp(−x/2) / [1 − exp(−x)]": r"z_{\mathrm{vib}}=\frac{\exp(-x/2)}{1-\exp(-x)}",
    "F_vib = RT Σ w_q {x/2 + ln[1 − exp(−x)]} / Z": r"F_{\mathrm{vib}}=\frac{RT}{Z}\sum_q w_q\left\{\frac{x}{2}+\ln\left[1-\exp(-x)\right]\right\}",
    "U_vib = RT Σ w_q x {1/2 + 1/[exp(x)−1]} / Z": r"U_{\mathrm{vib}}=\frac{RT}{Z}\sum_q w_qx\left\{\frac{1}{2}+\frac{1}{\exp(x)-1}\right\}",
    "S_vib = R Σ w_q {x/[exp(x)−1] − ln[1−exp(−x)]} / Z": r"S_{\mathrm{vib}}=\frac{R}{Z}\sum_q w_q\left\{\frac{x}{\exp(x)-1}-\ln\left[1-\exp(-x)\right]\right\}",
    "ΔF_screen = ΔE_static + ΔU_vib,proxy − TΔS_proxy": r"\Delta F_{\mathrm{screen}}=\Delta E_{\mathrm{static}}+\Delta U_{\mathrm{vib,proxy}}-T\Delta S_{\mathrm{proxy}}",
    "ΔF_screen = ΔE_static + ΔF_vib,proxy": r"\Delta F_{\mathrm{screen}}=\Delta E_{\mathrm{static}}+\Delta F_{\mathrm{vib,proxy}}",
    "N[f.u.] = N[simulation box] / (Z × n_a × n_b × n_c)": r"N[\mathrm{f.u.}]=\frac{N[\mathrm{simulation\ box}]}{Z\,n_a n_b n_c}",
    "μ(T,P) = μ_ref(T) + RT ln[f(T,P)/f_ref]": r"\mu(T,P)=\mu_{\mathrm{ref}}(T)+RT\ln\left[\frac{f(T,P)}{f_{\mathrm{ref}}}\right]",
    "f(T,P) = φ(T,P) P": r"f(T,P)=\phi(T,P)P",
    "Ω_ads,i(f,T) = −RT ∫(0→f) N_i(f′,T) dln f′": r"\Omega_{\mathrm{ads},i}(f,T)=-RT\int_{0}^{f}N_i(f',T)\,d\ln f'",
    "Ω_i = F_host,i(T) + P_ext V_i + Ω_ads,i": r"\Omega_i=F_{\mathrm{host},i}(T)+P_{\mathrm{ext}}V_i+\Omega_{\mathrm{ads},i}",
    "ΔΩ = ΔF_host(T) + P_ext ΔV − RT ∫(0→f) ΔN(f′,T) dln f′": r"\Delta\Omega=\Delta F_{\mathrm{host}}(T)+P_{\mathrm{ext}}\Delta V-RT\int_{0}^{f}\Delta N(f',T)\,d\ln f'",
    "PΔV[kJ/mol-f.u.] = 6.02214076×10⁻⁵ P[bar] ΔV[Å³/f.u.]": r"P\Delta V[\mathrm{kJ\,mol^{-1}\,f.u.^{-1}}]=6.02214076\times10^{-5}\,P[\mathrm{bar}]\,\Delta V[\mathrm{\mathring A^3\,f.u.^{-1}}]",
    "I(0→f_0) = ∫ K_H f dln f = K_H f_0": r"I(0\to f_0)=\int K_H f\,d\ln f=K_H f_0",
    "ΔΩ = ΔF_host − RT N_sat ln(1+bP)": r"\Delta\Omega=\Delta F_{\mathrm{host}}-RTN_{\mathrm{sat}}\ln(1+bP)",
    "P_eq = {exp[ΔF_host/(RTN_sat)] − 1}/b": r"P_{\mathrm{eq}}=\frac{\exp\left[\Delta F_{\mathrm{host}}/(RTN_{\mathrm{sat}})\right]-1}{b}",
    "Ω_ads,c = −RTK_s f": r"\Omega_{\mathrm{ads},c}=-RTK_s f",
    "ΔΩ_residual = ΔΩ_zero + RTK_s f": r"\Delta\Omega_{\mathrm{residual}}=\Delta\Omega_{\mathrm{zero}}+RTK_s f",
    "dΔΩ/dlnP ≈ −RT ΔN(P_eq)": r"\frac{d\Delta\Omega}{d\ln P}\approx-RT\,\Delta N(P_{\mathrm{eq}})",
    "δlnP_eq ≈ δ(ΔF_host) / [RT ΔN(P_eq)]": r"\delta\ln P_{\mathrm{eq}}\approx\frac{\delta(\Delta F_{\mathrm{host}})}{RT\,\Delta N(P_{\mathrm{eq}})}",
}


def escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def escape_url(url: str) -> str:
    return url.replace("%", r"\%").replace("#", r"\#").replace("_", r"\_")


def inline(text: str) -> str:
    text = text.replace("\u2013", "-").replace("\u2014", "-").replace("\u2011", "-")
    placeholders: list[str] = []

    def stash(value: str) -> str:
        token = f"@@LATEX{len(placeholders)}@@"
        placeholders.append(value)
        return token

    # Preserve explicit Markdown inline math as real TeX math.  The source
    # guide uses $...$ only for short symbols/equations; display equations are
    # still represented by the explicit [eq] mapping below.
    text = re.sub(r"\$([^$\n]+?)\$", lambda m: stash(rf"\({m.group(1)}\)"), text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
                  lambda m: stash(rf"\href{{{escape_url(m.group(2))}}}{{{escape(m.group(1))}}}"), text)
    text = re.sub(r"\*\*(.+?)\*\*", lambda m: stash(rf"\textbf{{{escape(m.group(1))}}}"), text)
    text = re.sub(r"`([^`]+)`", lambda m: stash(rf"\texttt{{{escape(m.group(1))}}}"), text)
    text = escape(text)
    for i, value in enumerate(placeholders):
        text = text.replace(f"@@LATEX{i}@@", value)
    return text


def markdown_to_tex(source: Path) -> str:
    lines = source.read_text(encoding="utf-8").splitlines()
    out = [r"\documentclass[UTF8,a4paper,11pt]{ctexart}",
           r"\usepackage{amsmath,amssymb,mathtools}",
           r"\usepackage{unicode-math}",
           r"\IfFontExistsTF{STIX Two Math}{\setmathfont{STIX Two Math}}{\setmathfont{Latin Modern Math}}",
           r"\IfFontExistsTF{Songti SC}{\setCJKmainfont{Songti SC}}{\setCJKmainfont{FandolSong}}",
           r"\IfFontExistsTF{PingFang SC}{\setCJKsansfont{PingFang SC}}{}",
           r"\usepackage[margin=2.0cm,top=2.25cm,bottom=1.75cm,headheight=15pt]{geometry}",
           r"\usepackage{xcolor,hyperref,fancyhdr,enumitem,booktabs}",
           r"\hypersetup{hidelinks,pdfcreator={XeLaTeX},pdftitle={Flexible MOF gate-opening fast screening principles}}",
           r"\urlstyle{same}", r"\setlength{\parindent}{2em}", r"\setlength{\parskip}{0.35em}",
           r"\linespread{1.06}", r"\emergencystretch=3em", r"\sloppy",
           r"\definecolor{navy}{HTML}{103B50}", r"\definecolor{teal}{HTML}{127C81}",
           r"\pagestyle{fancy}", r"\fancyhf{}", r"\fancyhead[L]{\small\sffamily FLEXIBLE-MOF / FAST SCREENING / PRINCIPLES}",
           r"\fancyhead[R]{\small\thepage}", r"\fancyfoot[C]{\small 条件性平衡预测 · 无 MD / DFT · 代理误差独立验证}",
           r"\renewcommand{\headrulewidth}{0.4pt}", r"\renewcommand{\footrulewidth}{0.2pt}",
           r"\makeatletter\renewcommand{\section}{\@startsection{section}{1}{0pt}{1.2ex plus .2ex}{.5ex}{\Large\bfseries\color{navy}}}\makeatother",
           r"\makeatletter\renewcommand{\subsection}{\@startsection{subsection}{2}{0pt}{1.0ex plus .2ex}{.35ex}{\large\bfseries\color{teal}}}\makeatother",
           r"\begin{document}"]
    in_items = False
    paragraph: list[str] = []

    def flush() -> None:
        nonlocal paragraph
        if paragraph:
            joined = paragraph[0]
            for part in paragraph[1:]:
                joined += (" " if joined[-1:].isascii() and part[:1].isascii() else "") + part
            out.append(inline(joined))
            paragraph = []

    def close_items() -> None:
        nonlocal in_items
        if in_items:
            out.append(r"\end{itemize}")
            in_items = False

    for raw in lines:
        line = raw.strip().replace("\u2013", "-").replace("\u2014", "-").replace("\u2011", "-")
        if not line:
            flush(); continue
        if line == "---PAGE---":
            flush(); close_items(); out.append(r"\clearpage"); continue
        if line.startswith("[eq] "):
            flush(); close_items()
            equation = FORMULAS.get(line[5:]) or FORMULAS.get(line[5:].replace("-", "–"))
            if equation is None:
                raise ValueError(f"unmapped equation: {line[5:]}")
            out.extend([r"\begin{equation*}", equation, r"\end{equation*}"])
            continue
        if line.startswith("# "):
            flush(); close_items(); out.append(rf"\section*{{{inline(line[2:])}}}"); continue
        if line.startswith("## "):
            flush(); close_items(); out.append(rf"\subsection*{{{inline(line[3:])}}}"); continue
        if line.startswith("- "):
            flush()
            if not in_items:
                out.append(r"\begin{itemize}[leftmargin=1.7em,itemsep=0.15em,topsep=0.15em]")
                in_items = True
            out.append(rf"\item {inline(line[2:])}")
            continue
        close_items(); paragraph.append(line)
    flush(); close_items(); out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def build(source: Path, out: Path, tex_out: Path | None = None, engine: str | None = None) -> Path:
    engine = engine or shutil.which("xelatex") or shutil.which("lualatex")
    if not engine:
        raise RuntimeError("XeLaTeX/LuaLaTeX is required to render the principles PDF")
    tex = markdown_to_tex(source)
    if tex_out is None:
        tex_out = out.with_suffix(".tex")
    tex_out.parent.mkdir(parents=True, exist_ok=True)
    tex_out.write_text(tex, encoding="utf-8")
    out.parent.mkdir(parents=True, exist_ok=True)
    build_dir = out.parent / ".latex-build"
    if build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir()
    local_tex = build_dir / tex_out.name
    local_tex.write_text(tex, encoding="utf-8")
    command = [engine, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error", local_tex.name]
    for _ in range(2):
        subprocess.run(command, cwd=build_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    rendered = build_dir / tex_out.with_suffix(".pdf").name
    if not rendered.is_file():
        raise RuntimeError("XeLaTeX completed without a PDF")
    shutil.copy2(rendered, out)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tex-out", type=Path)
    parser.add_argument("--engine")
    args = parser.parse_args()
    print(build(args.source, args.out, args.tex_out, args.engine))
