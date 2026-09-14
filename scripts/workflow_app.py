#!/usr/bin/env python3
"""Local-only Streamlit upload/run UI. Run with streamlit run, not python."""
import json
from pathlib import Path
import subprocess
import sys
import uuid
import streamlit as st

ROOT=Path(__file__).resolve().parents[1]
st.set_page_config(page_title="Flexible MOF workflow",layout="wide")
st.title("柔性 MOF：两相 CIF → 框架自由能 → 吸附 → 平衡压力")
st.caption("Local only · 不上传到第三方服务 · P_eq 不等于实验开门/关门压力")
method=st.selectbox("框架自由能方法 / Host free energy",["快速筛选：静态 ML 能量 + 弹性网络熵（不做 MD/DFT）","有限位移声子：较慢，虚频会阻断"])
st.info("通用 CIF 必须附 clean-cif 审核记录与审计 JSON。当前内置自动修复仅支持示例中的公开无氢 MIL-53(Al) 4×4×4 文件。未知无序/溶剂不会自动放行。")
left,right=st.columns(2)
closed=left.file_uploader("闭孔相 closed.cif",type="cif")
opened=right.file_uploader("开孔相 open.cif",type="cif")
spec=st.file_uploader("工作流配置 workflow.json",type="json")
support=st.file_uploader("配套文件：force_field.json、清洗记录、审计、气体模型等（扁平目录）",type=["json","cif"],accept_multiple_files=True)
st.markdown("配置格式及示例见 `references/end_to_end.md`。模型、PACMOF 和 RASPA 使用本机已安装路径；上传的 JSON 不是可执行脚本。")
if spec:
    try:
        config=json.loads(spec.getvalue());st.json(config,expanded=False)
    except Exception as e:
        st.error(str(e)); config=None
else: config=None
if st.button("保存输入并运行 / Import and run",disabled=not(closed and opened and config)):
    try:
        names=[f.name for f in support]
        if any(Path(n).name!=n or n in {"workflow.json","closed.cif","open.cif"} for n in names) or len(names)!=len(set(names)):
            raise ValueError("配套文件名不能重复、含路径或覆盖 workflow/open/closed")
        dest=ROOT/"output/ui_uploads"/uuid.uuid4().hex;dest.mkdir(parents=True)
        for name,file in (("closed.cif",closed),("open.cif",opened)):
            (dest/name).write_bytes(file.getvalue())
        for f in support:(dest/f.name).write_bytes(f.getvalue())
        config["phases"]["closed"]["cif"]="closed.cif"
        config["phases"]["open"]["cif"]="open.cif"
        config["host"]["method"]="elastic-network" if method.startswith("快速") else "phonons"
        (dest/"workflow.json").write_text(json.dumps(config,indent=2))
        log=(dest/"workflow.log").open("w")
        process=subprocess.Popen([sys.executable,str(ROOT/"scripts/run_workflow.py"),"--spec",str(dest/"workflow.json"),"--output-dir",str(dest/"run")],stdout=log,stderr=subprocess.STDOUT)
        log.close()
        st.session_state["job"]=(process,str(dest))
    except Exception as e:st.error(str(e))

@st.fragment(run_every=3)
def progress():
    if "job" not in st.session_state:return
    process,directory=st.session_state["job"];directory=Path(directory)
    st.write("结果目录 / Results:",str(directory))
    code=process.poll()
    st.write("运行中 / Running" if code is None else f"已退出 / Exited: {code} (2 = prediction gated; inspect report)")
    if (directory/"workflow.log").exists():st.code((directory/"workflow.log").read_text()[-16000:])
    report=directory/"run/workflow_summary.json"
    if report.exists():
        st.json(json.loads(report.read_text()),expanded=False)
        st.download_button("下载运行摘要 / Download summary",report.read_bytes(),"workflow_summary.json",mime="application/json")
progress()
