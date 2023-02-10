import logging

import pandas as pd

from datetime import date

from chebi_tools.ChEBIStandardizer import ChEBIStandardizer

import gradio as gr



STD = ChEBIStandardizer()


def today():
    return date.today().strftime("%y%m%d")


def process(names):
    return STD.process_many(names.split('\n'))


def export_csv(df):
    fn = f"{today()}-metabolomics-standardizer-output.csv"
    df.to_csv(fn, index=False)
    return gr.File.update(value=fn, visible=True)


interface =  gr.Blocks()
with interface:

    gr.Markdown("Metabolomics Standardizer")
    btn = gr.Button(value="Start")
    text_input = gr.Textbox(label="Input", max_lines=1000, interactive=True, lines=10, placeholder='Mevalonate\nmevalonic acid')
    export = gr.Button("Export")
    df_output = gr.Dataframe(label="Output", interactive=False)
    btn.click(process, inputs=text_input, outputs=df_output)

    csv = gr.File(interactive=False, visible=False)
    export.click(export_csv, df_output, csv)

if __name__ == '__main__':
    interface.launch()   