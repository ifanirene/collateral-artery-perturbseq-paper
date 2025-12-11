"""
/**
 * @description 
 * Generates an interactive, searchable HTML report for topic annotations with improved UI/UX.
 * It reads a summary CSV and per-topic markdown annotations, then renders a modern page with
 * search, sort, expand/collapse controls, a table of contents, badges, anchors, and copy-link.
 * 
 * Key features:
 * - Sticky toolbar with search, sort (Topic ID / Name), expand/collapse all
 * - Sidebar table of contents with quick jump links to topics
 * - Keyword and gene badges, anchor links per topic, copy-link button
 * - Responsive layout, refined typography, and subtle animations
 * 
 * @dependencies
 * - pandas: read summary CSV
 * - markdown: render markdown content into HTML
 * 
 * @examples
 * - python 04_generate_html_report.py summary.csv annotations_dir report.html
 */
"""

import pandas as pd
import argparse
import os
from markdown import markdown
from datetime import datetime
import json

def generate_html_overview(summary_csv, annotations_dir, output_html):
    """
    Generates a single HTML file to display topic summaries and annotations with a modern layout.
    """
    try:
        summary_df = pd.read_csv(summary_csv)
    except FileNotFoundError:
        print(f"Error: Summary CSV not found at {summary_csv}")
        return

    generated_on = datetime.now().strftime('%Y-%m-%d %H:%M')
    num_topics = len(summary_df)

    html_template = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Topic Annotation Overview</title>
        <style>
            :root {{
                /* Light theme (default) */
                --bg: #ffffff;
                --panel: #fafafa;
                --text: #111827;
                --muted: #6b7280;
                --accent: #2563eb;
                --accent-2: #059669;
                --border: #e5e7eb;
                --badge-bg: #f3f4f6;
            }}
            body.theme-dark {{
                --bg: #0b0c10;
                --panel: #121317;
                --text: #e9ecef;
                --muted: #a7b1c2;
                --accent: #3b82f6;
                --accent-2: #10b981;
                --border: #1f2430;
                --badge-bg: #1b2030;
            }}
            * {{ box-sizing: border-box; }}
            body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; line-height: 1.6; margin: 0; background: var(--bg); color: var(--text); }}
            .page {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
            .header {{ position: sticky; top: 0; z-index: 10; background: linear-gradient(180deg, rgba(11,12,16,0.95) 0%, rgba(11,12,16,0.75) 100%); backdrop-filter: blur(6px); border-bottom: 1px solid var(--border); }}
            .header-inner {{ max-width: 1400px; margin: 0 auto; padding: 12px 20px; display: grid; grid-template-columns: 1fr 2fr; gap: 12px; align-items: center; }}
            .title {{ font-size: 20px; font-weight: 700; letter-spacing: 0.3px; }}
            .stats {{ display: flex; gap: 10px; color: var(--muted); font-size: 12px; margin-top: 4px; }}
            .controls {{ display: grid; grid-template-columns: 1.5fr 0.9fr auto; gap: 10px; }}
            .input {{ width: 100%; padding: 10px 12px; border: 1px solid var(--border); background: var(--panel); color: var(--text); border-radius: 8px; outline: none; transition: border-color 0.2s; }}
            .input:focus {{ border-color: var(--accent); }}
            .select {{ padding: 10px 12px; border: 1px solid var(--border); background: var(--panel); color: var(--text); border-radius: 8px; }}
            .btn {{ padding: 10px 12px; background: var(--panel); border: 1px solid var(--border); color: var(--text); border-radius: 8px; cursor: pointer; transition: transform 0.05s ease, border-color 0.2s; }}
            .btn:hover {{ border-color: var(--accent); }}
            .btn:active {{ transform: translateY(1px); }}
            .layout {{ display: grid; grid-template-columns: 260px 1fr; gap: 16px; margin-top: 16px; }}
            .sidebar {{ position: sticky; top: 64px; height: calc(100vh - 80px); overflow: auto; padding: 12px; border: 1px solid var(--border); border-radius: 12px; background: var(--panel); }}
            .toc-title {{ font-size: 13px; color: var(--muted); margin-bottom: 6px; }}
            .toc a {{ display: block; color: var(--text); text-decoration: none; padding: 6px 8px; margin: 2px 0; border-radius: 6px; border: 1px solid transparent; font-size: 13px; }}
            .toc a:hover {{ background: var(--badge-bg); border-color: var(--border); }}
            .content {{ border: 1px solid var(--border); border-radius: 12px; background: var(--panel); padding: 16px; min-height: 70vh; }}
            .topic-title {{ margin-top: 0; font-size: 18px; letter-spacing: 0.2px; }}
            .mini-btn {{ padding: 6px 8px; background: var(--badge-bg); border: 1px solid var(--border); color: var(--text); border-radius: 6px; cursor: pointer; font-size: 12px; }}
            .topic-content h3 {{ margin-top: 16px; border-bottom: 1px solid var(--border); padding-bottom: 6px; }}
            .badges {{ display: flex; flex-wrap: wrap; gap: 6px; margin: 8px 0 12px; }}
            .badge {{ display: inline-block; background: var(--badge-bg); border: 1px solid var(--border); color: var(--text); padding: 4px 8px; font-size: 12px; border-radius: 999px; }}
            .badge.accent {{ border-color: var(--accent); color: #cfe1ff; }}
            .footer {{ text-align: center; color: var(--muted); font-size: 12px; margin: 16px 0 6px; }}
            .back-to-top {{ position: fixed; right: 16px; bottom: 16px; z-index: 20; }}
            @media (max-width: 1024px) {{ .layout {{ grid-template-columns: 1fr; }} .sidebar {{ display: none; }} .header-inner {{ grid-template-columns: 1fr; }} .controls {{ grid-template-columns: 1fr 1fr; }} }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="header-inner">
                <div>
                    <div class="title">Topic Annotation Overview</div>
                    <div class="stats">
                        <div>{{NUM_TOPICS}} topics</div>
                        <div>Generated: {{GENERATED_ON}}</div>
                    </div>
                </div>
                <div class="controls">
                    <input id="searchInput" class="input" type="text" onkeyup="filterTopics()" placeholder="Search by topic, keyword, gene, or text...">
                    <select id="sortSelect" class="select" onchange="sortTopics()">
                        <option value="topic">Sort: Topic ID</option>
                        <option value="name_asc">Sort: Name A→Z</option>
                        <option value="name_desc">Sort: Name Z→A</option>
                    </select>
                     <button class="btn" onclick="toggleTheme()" id="themeBtn">Dark mode</button>
                </div>
            </div>
        </div>
        <div class="page">
            <div class="layout">
                 <aside class="sidebar">
                     <div class="toc-title">Topics</div>
                     <div class="toc" id="toc">
                         {{TOC_HTML}}
                     </div>
                 </aside>
                 <main class="content">
                     <div id="topic-detail" class="topic-content">
                         <p style="color: var(--muted)">Select a topic on the left to view details.</p>
                     </div>
                 </main>
            </div>
            <div class="footer">End of report</div>
        </div>
        <button class="btn back-to-top" onclick="window.scrollTo({{top:0, behavior:'smooth'}})">Back to top</button>

        <script>
            function selectTopic(id) {{
                const data = window.TOPIC_DATA[id];
                const container = document.getElementById('topic-detail');
                if (!data) {{ container.innerHTML = '<p style="color: var(--muted)">No data.</p>'; return; }}

                // Build HTML for selected topic
                const figDir = '../enrichment_figures';
                const keggFig = `${figDir}/program_${id}_kegg_enrichment.png`;
                const procFig = `${figDir}/program_${id}_process_enrichment.png`;

                const figures = `
                    <div>
                        <h3>Enrichment Figures</h3>
                        <div style="display:grid;grid-template-columns:1fr 1fr;gap:12px;">
                            <div><img src="${keggFig}" alt="KEGG enrichment" style="max-width:100%;border:1px solid var(--border);border-radius:8px;"></div>
                            <div><img src="${procFig}" alt="Process enrichment" style="max-width:100%;border:1px solid var(--border);border-radius:8px;"></div>
                        </div>
                    </div>`;

                container.innerHTML = `
                    <h2 class="topic-title" id="topic-${id}">Topic ${id}: ${data.name}</h2>
                    <div class='keywords'><strong>Keywords:</strong> <div class='badges'>${data.keywords_html}</div></div>
                    <h3>Summary</h3>
                    <p>${data.summary}</p>
                    ${data.top_genes_block}
                    <h3>Full LLM Annotation</h3>
                    ${data.full_annotation}
                    ${figures}
                `;
            }

            function filterTopics() {{
                const input = document.getElementById('searchInput');
                const filter = input.value.toUpperCase();
                const toc = document.getElementById('toc');
                const links = toc.getElementsByTagName('a');
                for (let i = 0; i < links.length; i++) {{
                    const txtValue = links[i].textContent || links[i].innerText;
                    if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                        links[i].style.display = '';
                    }} else {{
                        links[i].style.display = 'none';
                    }}
                }}
            }}

            function toggleTheme() {{
                const body = document.body;
                const btn = document.getElementById('themeBtn');
                body.classList.toggle('theme-dark');
                const isDark = body.classList.contains('theme-dark');
                btn.innerText = isDark ? 'Light mode' : 'Dark mode';
                try {{ localStorage.setItem('topic_theme', isDark ? 'dark' : 'light'); }} catch (e) {{}}
            }}
            (function initTheme() {{
                try {{
                    const pref = localStorage.getItem('topic_theme');
                    if (pref === 'dark') {{ document.body.classList.add('theme-dark'); document.getElementById('themeBtn').innerText = 'Light mode'; }}
                }} catch (e) {{}}
            }})();

            function sortTopics() {{
                const select = document.getElementById('sortSelect');
                const value = select.value;
                const container = document.getElementById('toc');
                const items = Array.from(container.getElementsByTagName('a'));

                items.sort((a, b) => {{
                    if (value === 'topic') {{
                        const ida = parseInt(a.getAttribute('data-topic-id'));
                        const idb = parseInt(b.getAttribute('data-topic-id'));
                        return ida - idb;
                    }} else if (value === 'name_asc' || value === 'name_desc') {{
                        const na = a.getAttribute('data-topic-name').toLowerCase();
                        const nb = b.getAttribute('data-topic-name').toLowerCase();
                        const cmp = na.localeCompare(nb);
                        return value === 'name_asc' ? cmp : -cmp;
                    }}
                    return 0;
                }});

                items.forEach(it => container.appendChild(it));
            }}

            function copyLink(id) {{
                const url = window.location.origin + window.location.pathname + '#topic-'+id;
                navigator.clipboard.writeText(url).then(() => {{
                    const btn = document.getElementById('copy-topic-'+id);
                    if (btn) {{
                        const original = btn.innerText;
                        btn.innerText = 'Copied!';
                        setTimeout(()=> btn.innerText = original, 1200);
                    }}
                }});
            }}
        </script>
    </body>
    </html>
    """

    topics_html_content = ""
    toc_html_content = ""
    topic_data_js = []
    print("Generating HTML content for each topic...")
    for _, row in summary_df.iterrows():
        topic_id = row['Topic']
        topic_name = row['Name']
        keywords = row['Keywords']
        summary = row['Summary']

        # Top genes badges
        top_genes_block = ""
        if 'Top_Genes' in row and pd.notna(row['Top_Genes']):
            top_genes_list = [g.strip() for g in str(row['Top_Genes']).split(',') if g.strip()]
            top_gene_badges = "".join([f"<span class='badge accent'>{g}</span>" for g in top_genes_list[:10]])
            top_genes_block = f"""
                <div class='genes'>
                    <h3>Top Genes</h3>
                    <div class='badges'>{top_gene_badges}</div>
                </div>
            """

        # Read full annotation
        annotation_file = os.path.join(annotations_dir, f'topic_{topic_id}_annotation.md')
        full_annotation = ""
        if os.path.exists(annotation_file):
            with open(annotation_file, 'r', encoding='utf-8') as f:
                full_annotation = markdown(f.read())
        else:
            print(f"  ✗ Annotation file not found for Topic {topic_id}")

        # TOC entry (click to load)
        toc_html_content += f"<a href='javascript:void(0)' onclick='selectTopic({topic_id})' data-topic-id='{topic_id}' data-topic-name='{topic_name}'>Topic {topic_id}: {topic_name}</a>"

        # Keyword badges
        keyword_badges = ""
        if isinstance(keywords, str) and keywords.strip():
            keyword_list = [k.strip() for k in keywords.split(',') if k.strip()]
            keyword_badges = "".join([f"<span class='badge'>{k}</span>" for k in keyword_list])

        # Prepare per-topic data for JS single-view rendering
        topic_data_js.append({
            'id': int(topic_id),
            'name': topic_name,
            'keywords_html': keyword_badges,
            'summary': summary,
            'top_genes_block': top_genes_block,
            'full_annotation': full_annotation,
        })
        print(f"  ✓ Processed Topic {topic_id}")

    final_html = html_template
    final_html = final_html.replace("{{TOC_HTML}}", toc_html_content)
    final_html = final_html.replace("{{NUM_TOPICS}}", str(num_topics))
    final_html = final_html.replace("{{GENERATED_ON}}", generated_on)
    # Normalize any double-brace artifacts back to single braces so CSS/JS render properly
    final_html = final_html.replace("{{", "{").replace("}}", "}")

    # Inject per-topic data for client-side rendering using JSON to ensure valid JS strings
    topic_data_script = "<script>\nwindow.TOPIC_DATA = {};\n"
    for d in topic_data_js:
        payload = json.dumps({
            'name': d['name'],
            'keywords_html': d['keywords_html'],
            'summary': d['summary'],
            'top_genes_block': d['top_genes_block'],
            'full_annotation': d['full_annotation'],
        })
        topic_data_script += f"window.TOPIC_DATA[{d['id']}] = {payload};\n"
    topic_data_script += "</script>"

    final_html = final_html.replace("</body>", topic_data_script + "\n</body>")
    try:
        with open(output_html, 'w', encoding='utf-8') as f:
            f.write(final_html)
        print(f"\n✓ Successfully generated HTML overview at: {output_html}")
    except IOError as e:
        print(f"  ✗ Error writing HTML file: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate an HTML overview from topic summaries and annotations.")
    parser.add_argument("summary_csv", help="Path to the topic summary CSV file.")
    parser.add_argument("annotations_dir", help="Directory containing individual topic markdown annotation files.")
    parser.add_argument("output_html", help="Path to the output HTML file.")
    args = parser.parse_args()

    generate_html_overview(args.summary_csv, args.annotations_dir, args.output_html) 