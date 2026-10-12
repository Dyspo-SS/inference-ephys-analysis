"""Command-line entry points; all paths are explicit and no results are hard coded."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback
from datetime import datetime,timezone
import pandas as pd
from . import __version__
from . import io,rates,rsa,glm,behavior,plotting
from .common import checksum,write_json

def manifest(config_path,c):
    versions={}
    for package in ["numpy","pandas","scipy","matplotlib"]:
        try:versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:versions[package]="unavailable"
    try:
        repo=Path(__file__).resolve().parents[2]
        commit=subprocess.run(["git","rev-parse","HEAD"],cwd=repo,capture_output=True,text=True,check=True).stdout.strip()
    except (OSError,subprocess.CalledProcessError):commit=None
    return dict(Status="running",Started_UTC=datetime.now(timezone.utc).isoformat(),Package_Version=__version__,Git_Commit=commit,
                Python=sys.version,Platform=platform.platform(),Dependencies=versions,Config_SHA256=checksum(config_path),Configuration=c)

def output_inventory(out):
    rows=[]
    for path in sorted(out.rglob("*")):
        if not path.is_file() or path.name in ["output_inventory.json","run_manifest.json"]:continue
        row={"File":path.relative_to(out).as_posix(),"Bytes":path.stat().st_size,"SHA256":checksum(path)}
        if path.suffix==".csv":row["Columns"]=list(pd.read_csv(path,nrows=0).columns)
        rows.append(row)
    write_json(rows,out/"audit/output_inventory.json")

def main(argv=None):
    parser=argparse.ArgumentParser(description="Animal-level electrophysiology analysis")
    parser.add_argument("command",choices=["validate","run","plot"])
    parser.add_argument("--config",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--no-plots",action="store_true",help="Export analytical outputs without figures")
    args=parser.parse_args(argv);c,cp=io.read_config(args.config);out=args.out.resolve()
    if args.command=="plot":
        record=out/"audit/run_manifest.json"
        if not record.exists():parser.error("plot requires an existing analyzed output directory")
        previous=json.loads(record.read_text(encoding="utf-8"))
        if previous.get("Status")!="complete":parser.error("The analysis did not complete")
        if previous.get("Config_SHA256")!=checksum(cp):parser.error("Configuration changed; rerun analysis into a new directory")
        plotting.run(out,c);output_inventory(out);return
    if out.exists() and any(out.iterdir()):parser.error("Output directory must be empty to prevent mixing analysis runs")
    out.mkdir(parents=True,exist_ok=True)
    for folder in ["audit","tables","statistics","plotdata","figures"]:(out/folder).mkdir(exist_ok=True)
    rec=manifest(cp,c);write_json(rec,out/"audit/run_manifest.json");write_json(c,out/"audit/config_used.json")
    try:
        d,source=io.load_counts(c,cp,out);rec["Input_SHA256"]=checksum(source)
        if args.command=="validate":rec["Status"]="validated_only"
        else:
            u,_,ar=rates.run(d,c,out)
            if c["glm"]["enabled"]:glm.run(d,c,out)
            if c["rsa"]["enabled"]:
                rsa.run(d,c,out,u)
                if c["study"]=="intervention":rsa.count_controls(d,c,out,ar)
            behavior.run(c,cp,out)
            if c.get("behavior",{}).get("input"):
                rec["Behavior_Input_SHA256"]=checksum((cp.parent/c["behavior"]["input"]).resolve())
            if not args.no_plots:plotting.run(out,c)
            rec["Status"]="complete"
        rec["Completed_UTC"]=datetime.now(timezone.utc).isoformat();write_json(rec,out/"audit/run_manifest.json");output_inventory(out)
    except Exception as exc:
        rec["Status"]="failed";rec["Error"]=str(exc);write_json(rec,out/"audit/run_manifest.json")
        (out/"audit/error.txt").write_text(traceback.format_exc(),encoding="utf-8")
        raise
    print(f"{rec['Status']}: {out}",flush=True)

if __name__=="__main__":main()
