//Draw the secondary protein structures
function drawStructurePan(structure) {
    const canvas = document.getElementById('proteinStructure-pan');
    const ctx = canvas.getContext('2d');

    const legendWidth = 200;
    const width = (canvas.width - legendWidth) / alignment_length;

    // Draw the structure bars
    for (let i = 0; i < structure.length; i++) {
        let color;
        switch (structure[i]) {
            case 'H':
                color = 'red';
                break;
            case 'E':
                color = 'green';
                break;
            default:
                color = 'white';
        }
            ctx.fillStyle = color;
            //Add 1 for 0 offset, then substract 1 to start at origin
            ctx.fillRect((X2_to_X_array[i + 1] - 1) * width, 0, width, canvas.height);
    }

    // Draw the legend on the far right
    const legendData = [
        {color: 'red', label: 'Helix'},
        {color: 'green', label: 'Sheet'},
    ];

    const legendBoxSize = 30; // Size of each colored box in the legend
    const padding = 10; // Space between each legend item
    let yOffset = 5; // Initial vertical position for the first legend item
    let xOffset = 20;

    ctx.font = '15px Arial';

    for (const item of legendData) {
        ctx.fillStyle = item.color;
        ctx.fillRect(canvas.width - legendWidth + padding + xOffset, yOffset, legendBoxSize, legendBoxSize);

        ctx.fillStyle = 'black'; // Text color for the labels
        ctx.fillText(item.label, canvas.width - legendWidth + 2 * padding + legendBoxSize + xOffset, yOffset + legendBoxSize - 5);

        yOffset += legendBoxSize + padding ;
    }
}

function renderHeatmapPan(data) {
    const heatmapEl = document.getElementById('heatmap-pan');
    const cellWidth = window_length / alignment_length;
    const cellHeight = 5;
    const fullHeatmapEl = document.getElementById('full-heatmap-pan');
    const container = document.createElement('div');
    const fragment = document.createDocumentFragment(); // Create a document fragment

    data.forEach(cellData => {
        const xPosition = +cellData["X"];  // + ensures the string is converted to a number
        const yPosition = +cellData["Y"];
        //const score = +cellData["Score"];

        const cell = document.createElement('div');
        cell.style.width = (cellWidth + 0.5) + 'px';
        cell.style.height = cellHeight + 'px';

        const wt = String(cellData["WT"] ?? "").trim();
        const sub = String(cellData["Sub"] ?? "").trim();
        const scoreRaw = cellData["Score"];
        const score = (scoreRaw === "" || scoreRaw == null) ? NaN : +scoreRaw;

        const isGap = (wt.startsWith("-") || sub.startsWith("-") || scoreRaw === "" || scoreRaw == null || Number.isNaN(score));
        cell.style.backgroundColor = isGap ? "#e6e6e6" : colorScalePan(score);

        // (Recommended) full heatmap should usually NOT render letters for performance.
        // If you want letters only in zoom, leave this blank:
        cell.innerText = "";


        cell.style.position = 'absolute';
        cell.style.left = cellWidth * (xPosition - 1) + 'px';  // -1 because X is 1-indexed
        cell.style.top = cellHeight * (yPosition - 1) + 'px';  // -1 because Y is 1-indexed
        cell.dataset.x = cellData["X"];

        const tooltip_h = document.createElement("div");
            tooltip_h.classList.add("tooltip");
            document.body.appendChild(tooltip_h);

        cell.addEventListener("mouseover", (e) => {
           tooltip_h.style.visibility = "visible";
           tooltip_h.innerHTML = `
               B73 Position: ${cellData["X2"]}<br>
               Target Position: ${cellData["X3"]}<br>
               Genome: ${GN_array[cellData["Y"]]}<br>
               G.M.: ${GM_array[cellData["Y"]]}<br>
               Substitution: ${cellData["WT"]} to ${cellData["Sub"]}<br>
               Score: ${cellData["Score"]}
           `;
           tooltip_h.style.left = (e.pageX + 10) + "px";
           tooltip_h.style.top = (e.pageY + 10) + "px";
       });

       // Mouse move event (so the tooltip_h follows the mouse)
        cell.addEventListener("mousemove", (e) => {
            tooltip_h.style.left = (e.pageX + 10) + "px";
            tooltip_h.style.top = (e.pageY + 10) + "px";
        });

        // Mouse out event
        cell.addEventListener("mouseout", () => {
            tooltip_h.style.visibility = "hidden";
        });

        fragment.appendChild(cell); // Append to the fragment instead of directly to the DOM
    });

    fullHeatmapEl.appendChild(fragment);

    const colors = [
        '#00429d', '#3860aa', '#587fb3', '#78a0b7', '#9ac0b3', '#c1e19e',
        '#ffff00', '#ffd337', '#fea447', '#f1784d', '#db4c4d', '#bd2147', '#93003a'
    ];

    container.style.position = 'absolute';
    container.style.left = (window_length + 210) + 'px';

    const benignEffectLabel = document.createElement('div');
    benignEffectLabel.innerText = 'Benign effect';
    benignEffectLabel.className = 'labelCB';
    benignEffectLabel.style.width = '100px';
    benignEffectLabel.style.fontFamily = 'Arial';
    benignEffectLabel.style.fontSize = '15px';

    container.appendChild(benignEffectLabel);

    colors.forEach((color, index) => {
        const row = document.createElement('div');
        row.className = 'rowCB';

        const box = document.createElement('div');
        box.className = 'colorBox';
        box.style.backgroundColor = color;
        row.appendChild(box);

        const value = document.createElement('div');
        value.className = 'valueCB';
        value.style.fontFamily = 'Arial';
        value.style.fontSize = '15px';

        if(index == 0)
        {
            value.innerText = '> 0';
        } else if(index == 12)
        {
            value.innerText = '< -22';
        } else {
            value.innerText = (2 - (index * 2)) + '';
        }

        row.appendChild(value);

        if(index == 5)
        {
            const mildEffectLabel = document.createElement('div');
            mildEffectLabel.innerText = 'Mild effect';
            mildEffectLabel.className = 'labelCB';
            mildEffectLabel.style.width = '100px';
            mildEffectLabel.style.fontFamily = 'Arial';
            mildEffectLabel.style.fontSize = '15px';
            container.appendChild(mildEffectLabel);
        }

        container.appendChild(row);
    });

    // Strong Effect label
    const strongEffectLabel = document.createElement('div');
    strongEffectLabel.innerText = 'Strong Effect';
    strongEffectLabel.className = 'labelCB';
    strongEffectLabel.style.width = '100px';
    strongEffectLabel.style.fontFamily = 'Arial';
    strongEffectLabel.style.fontSize = '15px';
    strongEffectLabel.style.marginBottom = '8px';   // slight margin under label
    container.appendChild(strongEffectLabel);

    // Row wrapper (controls spacing to the next legend block)
    const rowL = document.createElement('div');
    rowL.className = 'rowCB';
    rowL.style.marginBottom = '30px';               // 30px before next feature

    // Grey box (make it a flex container so text centers)
    const boxL = document.createElement('div');
    boxL.className = 'colorBox';
    boxL.style.backgroundColor = '#DDDDDD';
    boxL.style.width = '80px';                      // adjust to taste
    boxL.style.height = '22px';                     // adjust to taste
    boxL.style.display = 'flex';
    boxL.style.alignItems = 'center';
    boxL.style.justifyContent = 'center';
    boxL.style.fontFamily = 'Arial';
    boxL.style.fontSize = '15px';
    boxL.style.color = '#000';
    boxL.innerText = 'INDEL';                       // text inside the box

    rowL.appendChild(boxL);
    container.appendChild(rowL);

    const containerWrapper = document.createElement('div');
    containerWrapper.style.border = '1px solid black';
    containerWrapper.style.borderRadius = '5px'; // Add rounded corners with a 5px radius
    //containerWrapper.style.margin = '15px'; // Add a 5px margin

    // Create a wrapper div
    const container2= document.createElement('div');
    container2.style.margin = '15px'; // Add a 5px margin to the wrapper

    ///This is custum code to add a legend for the color-coding of the genomes

    let GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Heterotic group';
    GenomeLabel.className = 'labelCB';
    GenomeLabel.style.width = '125px';
    GenomeLabel.style.fontFamily = 'Arial';
    GenomeLabel.style.fontSize = '15px';
    GenomeLabel.style.marginTop = '0px'; // Add margin here
    GenomeLabel.style.textDecoration = 'underline';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Stiff stalk';
    GenomeLabel.style.color = 'black';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Mix';
    GenomeLabel.style.color = '#666666';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Non-stiff-stalk';
    GenomeLabel.style.color = '#455edd';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Iodent';
    GenomeLabel.style.color = '#a807ed';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Lancaster';
    GenomeLabel.style.color = '#da9af5';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'European flint';
    GenomeLabel.style.color = '#9ae6f5';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Chinese';
    GenomeLabel.style.color = '#f50707';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Tang SiPingTou';
    GenomeLabel.style.color = '#fa7d7d';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Popcorn';
    GenomeLabel.style.color = '#ce58ce';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Sweet corn';
    GenomeLabel.style.color = 'pink';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Tropical';
    GenomeLabel.style.color = '#30c727';
    container2.appendChild(GenomeLabel);

    //Add this code back in once the PanAnd genomes are available
    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'PanAnd';
    GenomeLabel.style.color = '#773510';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Teosinte';
    GenomeLabel.style.color = '#ca854c';
    container2.appendChild(GenomeLabel);

    GenomeLabel = document.createElement('div');
    GenomeLabel.innerText = 'Hi/Lo';
    GenomeLabel.style.color = '#1B9E77';
    container2.appendChild(GenomeLabel);

    containerWrapper.appendChild(container2);

    container.appendChild(containerWrapper);

    fullHeatmapEl.appendChild(container);
}

function updateHighlightBoxPan(start) {
    const highlightBox = document.getElementById('highlight-box-pan');
    const zoomedWidth = 50;  // We're zooming into 50 units
    highlightBox.style.left = (start - 1) * scaleFactorPan + 'px';  // -1 because we're 1-indexed
    highlightBox.style.width = (zoomedWidth * scaleFactorPan - 3) + 'px';   //adjust for border length (3px)
}


function renderDomainsPan(geneModelElement, domainsArray) {
    // Clear old render + remove old tooltips
    geneModelElement.innerHTML = "";
    document.querySelectorAll(".tooltipdomain").forEach(t => t.remove());
    currentlyVisibleTooltip = null;

    const domainHeight = 20;
    const domainMarginTop = 5;

    if (domainsArray.length > 0) {
        geneModelElement.style.height = domainHeight + (domainHeight + domainMarginTop) * (domainsArray.length - 1) + 'px';
    } else {
        geneModelElement.style.height = '20px';
    }

    // Adding a light grey horizontal line
     const line = document.createElement('div');
     line.style.position = 'absolute';
     line.style.top = '0px'; // adjust this value based on your needs
     line.style.left = '0';
     line.style.width = '100%';
     line.style.height = '2px'; // height of the line
     line.style.backgroundColor = 'lightgrey';

     geneModelElement.appendChild(line);

    let currentTopPosition = 0;
    let PGM_height = 100;
    let colorSave = [];

    domainsArray.forEach((domain, index) => {
        const div = document.createElement('div');
        const tooltip = document.createElement("div");
        tooltip.classList.add("tooltipdomain");
        tooltip.style.zIndex = "1000"; // Ensure tooltip is on top
        tooltip.style.pointerEvents = "auto"; // Make sure you can interact with the tooltip
        document.body.appendChild(tooltip);

        div.className = 'domain';
        div.style.left = scalePan(domain.start2) + 'px';
        div.style.width = (scalePan(domain.end2) - scalePan(domain.start2)) + 'px';
        div.style.top = (domainHeight + domainMarginTop) * 2 * index + 'px';
        div.style.backgroundColor = pastelColors[index % pastelColors.length]; // Color assignment

        if (domain.pfam_id) {

                if(colorSave[domain.pfam_id])
                {
                    div.style.backgroundColor = colorSave[domain.pfam_id];
                } else {
                    colorSave[domain.pfam_id] = pastelColors[index % pastelColors.length];
                }
        }

        if (domain.pfam_name) {
            div.addEventListener("click", (e) => {
                // Check if the tooltip being clicked is the same as the currently visible one
                if (currentlyVisibleTooltip && currentlyVisibleTooltip === tooltip) {
                    // Toggle visibility
                    if (tooltip.style.visibility === "visible") {
                        tooltip.style.visibility = "hidden";
                        currentlyVisibleTooltip = null;
                    } else {
                        tooltip.style.visibility = "visible";
                        tooltip.innerHTML = `
                            Pfam id: <a target="_blank" href="https://www.ebi.ac.uk/interpro/entry/pfam/${domain.pfam_id}">${domain.pfam_id}</a><br>
                            Pfam name: ${domain.pfam_name}<br>
                            Interprot id: ${domain.inter_id}<br>
                            Interprot name: ${domain.inter_name}<br>
                            GO terms: ${domain.go}<br>
                            Start: ${domain.start}<br>
                            End: ${domain.end}<br>
                        `;
                        tooltip.style.left = (e.pageX + 10) + "px";
                        tooltip.style.top = (e.pageY + 10) + "px";
                        currentlyVisibleTooltip = tooltip;
                    }
                } else {
                    // Hide the previously visible tooltip, if any
                    if (currentlyVisibleTooltip) {
                        currentlyVisibleTooltip.style.visibility = "hidden";
                    }
                    // Show the tooltip for the clicked div
                    tooltip.style.visibility = "visible";
                    tooltip.innerHTML = `
                        Pfam id: <a target="_blank" href="https://www.ebi.ac.uk/interpro/entry/pfam/${domain.pfam_id}">${domain.pfam_id}</a><br>
                        Pfam name: ${domain.pfam_name}<br>
                        Interprot id: ${domain.inter_id}<br>
                        Interprot name: ${domain.inter_name}<br>
                        GO terms: ${domain.go}<br>
                        Start: ${domain.start}<br>
                        End: ${domain.end}<br>
                    `;
                    tooltip.style.left = (e.pageX + 10) + "px";
                    tooltip.style.top = (e.pageY + 10) + "px";
                    currentlyVisibleTooltip = tooltip;
                }
            });
        }
        geneModelElement.appendChild(div);

        // Create a new div for the Pfam name
        const pfamNameDiv = document.createElement('div');
        pfamNameDiv.innerText = domain.pfam_name;
        pfamNameDiv.className = 'truncate';
        pfamNameDiv.style.left = scalePan(domain.start2) + 'px';
        pfamNameDiv.style.width = (scalePan(domain.end2) - scalePan(domain.start2)) + 'px';
        pfamNameDiv.style.top = (((domainHeight + domainMarginTop) * 2 * index) + 30) + 'px';

        PGM_height = (((domainHeight + domainMarginTop) * 2 * index) + 30);

        geneModelElement.appendChild(pfamNameDiv);  // Add Pfam name div as a child of the domain div

    });
        //adjust the height of the pfam conatainer
        document.getElementById('pfamGeneModel-pan').style.height = PGM_height + "px";
    }

    function loadDSSPPan() {
        let dssp_file = './dssp/' + protein + '.tsv';

        fetch(dssp_file)
        .then(response => {
            if (!response.ok) {
                throw new Error('DSSP file not found'); // Throw an error if response is not ok
            }

            return response.text();
        })
        .then(data => {
            if(data.includes("<html") || data.includes("<script"))
            {
                return false;
            }

            const structure = parseDSSP(data);
            drawStructurePan(structure);  // Call the function here with the parsed data
        })
        .catch(error => {
            document.getElementById("dssp-pan").innerHTML = '<div class="text-content"><br>There is no secondary structure information for this protein.</div>'
            console.log('Error fetching DSSP file:', error);
        });
    }

    // Function to create the number line for pan-view, offset for Indels
    function createNumberLinePan(container) {
      // Clear old ticks/labels
      container.innerHTML = "";

      // Tick/label at alignment column 1 (may map to '-' if query starts with a gap)
      const tick0 = document.createElement('div');
      tick0.className = 'tick';
      tick0.style.left = scalePan(1) + 'px';

      const label0 = document.createElement('div');
      label0.className = 'label';
      label0.style.left = scalePan(1) + 'px';

      const v0 = X_to_X2["1"];
      label0.innerText = (v0 === undefined || v0 === null || v0 === "") ? "-" : v0;

      container.appendChild(tick0);
      container.appendChild(label0);

      // Major ticks every 100 alignment columns
      for (let i = 99; i <= alignment_length - 1; i += 100) {
        const tick = document.createElement('div');
        tick.className = 'tick';
        tick.style.left = scalePan(i) + 'px';

        const label = document.createElement('div');
        label.className = 'label';
        label.style.left = scalePan(i) + 'px';

        const alnX = String(i + 1);  // alignment column (1-based)
        const v = X_to_X2[alnX];
        label.innerText = (v === undefined || v === null || v === "") ? "-" : v;

        container.appendChild(tick);
        container.appendChild(label);
      }

      // Last tick/label at last alignment column
      const tickLast = document.createElement('div');
      tickLast.className = 'tick-last';
      tickLast.style.left = scalePan(alignment_length - 1) + 'px';

      const labelLast = document.createElement('div');
      labelLast.className = 'label-last';
      labelLast.style.left = (scalePan(alignment_length - 1) + 20) + 'px';

      const vLast = X_to_X2[String(alignment_length)];
      labelLast.innerText = (vLast === undefined || vLast === null || vLast === "") ? "-" : vLast;

      container.appendChild(tickLast);
      container.appendChild(labelLast);
    }


    // Function to create the number line
    function createZoomedNumberLinePan(container,start) {
        // start is alignment column (1-based) for the current 50-col window
        container.innerHTML = "";

        for (let i = 1; i <= 49; i += 2) {
            const tick = document.createElement('div');
            tick.className = 'tick';
            tick.style.left = scaleZoomPan(i) + 'px';

            const label = document.createElement('div');
            label.className = 'label';
            label.style.left = scaleZoomPan(i) + 'px';

            const alnX = start + i;                 // alignment column
            const v = X_to_X2[String(alnX)];
            label.innerText = (v === undefined || v === null || v === "") ? "-" : v;

            container.appendChild(tick);
            container.appendChild(label);
        }
    }

    function updateHeatmapZoomPan(data, start) {
      const end = start + 49;
      const zoomedHeatmapEl = document.getElementById('zoomed-heatmap-pan');
      const cellWidth = window_length / 50;
      const cellHeight = 20;
      const fragment = document.createDocumentFragment();

      zoomedHeatmapEl.innerHTML = ''; // clear

      const ROW_OFFSET = 1; // ignore heatmap row 1 (B73_v5)

      // --- heatmap cells ---
      data.forEach(cellData => {
        const x = +cellData["X"];      // Option B: X is alignment column (1-based)
        if (x < start || x > end) return;

        //const y = +cellData["Y"];
        const yRaw = +cellData["Y"];
        const wt = String(cellData["WT"] ?? "").trim();
        const sub = String(cellData["Sub"] ?? "").trim();

        const scoreRaw = cellData["Score"];
        const score = (scoreRaw === "" || scoreRaw == null) ? NaN : +scoreRaw;
        const isGap = (wt.startsWith("-") || sub.startsWith("-") || scoreRaw === "" || scoreRaw == null || Number.isNaN(score));

        const cell = document.createElement('div');
        cell.className = 'heatmap-cell';
        cell.style.width = (cellWidth + 0.5) + 'px';
        cell.style.height = cellHeight + 'px';
        cell.style.position = 'absolute';
        cell.style.left = ((x - start) * cellWidth) + 'px';
        //cell.style.top = (cellHeight * (y - 1)) + 'px';

        const y = yRaw - ROW_OFFSET;     // shift up by 1
        if (y < 1) return;               // drop the first row entirely
        cell.style.top = (cellHeight * (y - 1)) + "px";

        cell.style.backgroundColor = isGap ? "#e6e6e6" : colorScalePan(score);

        // centered letter in the zoom view
        cell.innerText = sub;
        cell.style.color = "black";
        cell.style.textAlign = "center";
        cell.style.lineHeight = cellHeight + "px";
        cell.style.fontFamily = "monospace";
        cell.style.paddingLeft = "0px";

        // datasets for tooltip
        cell.dataset.X = cellData["X"];
        cell.dataset.X2 = cellData["X2"];
        cell.dataset.Y = cellData["Y"];
        cell.dataset.X3 = cellData["X3"];
        cell.dataset.WT = wt;
        cell.dataset.Sub = sub;
        cell.dataset.Score = cellData["Score"];

        fragment.appendChild(cell);
      });

      zoomedHeatmapEl.appendChild(fragment);

      // --- dropdown selector ---
      const cell0 = document.createElement('div');
      cell0.style.width = (cellWidth + 50) + 'px';
      cell0.style.height = cellHeight + 'px';
      cell0.style.backgroundColor = "transparent";
      cell0.style.position = 'absolute';
      cell0.style.left = (50 * cellWidth + 6) + 'px';
      cell0.style.top = '-38px';
      cell0.style.paddingLeft = "6px";
      cell0.style.zIndex = "20";

      if (view_type == "Species") {
        cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option selected value="Species">Species</option><option value="Protein">Gene Model</option></select>';
      } else {
        cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option value="Species">Species</option><option selected value="Protein">Gene Model</option></select>';
      }
      zoomedHeatmapEl.appendChild(cell0);

      const LABEL_ROW_OFFSET = 0; // ignore heatmap row 1 (B73_v5)
      let save_top_px = 0;

      for (let yRaw = 1 + LABEL_ROW_OFFSET; yRaw < GN_size; yRaw++) {
        const y = yRaw - LABEL_ROW_OFFSET; // display row index (1..)

        const cell = document.createElement("div");
        cell.style.position = "absolute";
        cell.style.left = (50 * cellWidth + 6) + "px";
        cell.style.top  = (cellHeight * (y - 1)) + "px";
        cell.style.height = cellHeight + "px";
        cell.style.width  = (cellWidth + 15) + "px";
        cell.style.backgroundColor = "transparent";
        cell.style.paddingLeft = "6px";
        cell.style.whiteSpace = "nowrap";

        // Use arrays with yRaw (because arrays are aligned to original Y indexing)
        let GN_name_text = String(GN_array[yRaw]);
        if (GM_array[yRaw].startsWith("Zm00014ba")) GN_name_text = "Mo17_v2";
        else if (GM_array[yRaw].startsWith("Zm00014a")) GN_name_text = "Mo17_v1";

        if (view_type === "Species") {
          cell.innerHTML =
            "<span class='species_option'>" + GN_name_text + "</span>" +
            "<span class='protein_option' style='display:none;'>" + String(GM_array[yRaw]) + "</span>";
        } else {
          cell.innerHTML =
            "<span class='species_option' style='display:none;'>" + GN_name_text + "</span>" +
            "<span class='protein_option' style='display:inline;'>" + String(GM_array[yRaw]) + "</span>";
        }

        cell.style.color = colorGenome(String(GN_array[yRaw]));
        zoomedHeatmapEl.appendChild(cell);

        save_top_px = cellHeight * y;  // first free row below last visible row
      }

      const save_top = save_top_px + "px";


      // --- bottom WT row (alignment space, gapped) ---
      for (let idx = 0; idx < 50; idx++) {
        const alnX = start + idx; // alignment column
        const aa = (X_to_WT && X_to_WT[String(alnX)]) ? X_to_WT[String(alnX)] : "-";

        const label = document.createElement('div');
        label.className = 'label_bottom';
        label.style.position = 'absolute';
        label.style.left = (idx * cellWidth) + 'px';
        label.style.top = save_top;
        label.style.width = cellWidth + 'px';
        label.style.height = cellHeight + 'px';
        label.style.textAlign = 'center';
        label.style.lineHeight = cellHeight + 'px';
        label.style.fontFamily = 'monospace';
        label.innerText = aa;

        zoomedHeatmapEl.appendChild(label);
      }
    }


    function updateHeatmapZoomPan_SAVE(data, start) {
  const end = start + 49;
  const zoomedHeatmapEl = document.getElementById('zoomed-heatmap-pan');
  const cellWidth = window_length / 50;
  const cellHeight = 20;
  const fragment = document.createDocumentFragment();

  zoomedHeatmapEl.innerHTML = ''; // Clear previous zoomed heatmap content

  // -----------------------------
  // Draw heatmap cells
  // -----------------------------
  data.forEach(cellData => {
    const x = +cellData["X"];

    if (x >= start && x <= end) {
      const cell = document.createElement('div');
      cell.className = 'heatmap-cell';

      const yPosition = +cellData["Y"];
      //const score = +cellData["Score"];

      // IMPORTANT: WT/Sub are letters; do NOT use unary +
      const wild_type = String(cellData["WT"] ?? "");
      const mutation  = String(cellData["Sub"] ?? "");

      cell.style.width = (cellWidth + 0.5) + 'px';
      cell.style.height = cellHeight + 'px';
      //cell.style.backgroundColor = colorScalePan(score);
      const wt = String(cellData["WT"] ?? "").trim();
      const sub = String(cellData["Sub"] ?? "").trim();
      const scoreRaw = cellData["Score"];
      const score = (scoreRaw === "" || scoreRaw == null) ? NaN : +scoreRaw;

      const isGap = (wt.startsWith("-") || sub.startsWith("-") || scoreRaw === "" || scoreRaw == null || Number.isNaN(score));

      cell.style.backgroundColor = isGap ? "#e6e6e6" : colorScalePan(score);

      cell.style.position = 'absolute';
      cell.style.left = ((x - start) * cellWidth) + 'px';
      cell.style.top  = (cellHeight * (yPosition - 1)) + 'px';

      // datasets
      cell.dataset.X = cellData["X"];
      cell.dataset.X2 = cellData["X2"];
      cell.dataset.Y = cellData["Y"];
      cell.dataset.X3 = cellData["X3"];
      cell.dataset.WT = wild_type;
      cell.dataset.Sub = mutation;
      cell.dataset.Score = cellData["Score"];

      cell.innerText = mutation;
      cell.style.paddingLeft = "6px";

      fragment.appendChild(cell);
    }
  });

  zoomedHeatmapEl.appendChild(fragment);

  // -----------------------------
  // Dropdown menu
  // -----------------------------
  const cell0 = document.createElement('div');
  let x = 50;                    // declare locally (avoid globals)
  let yPosition = 1;             // not used for cell0, but keep declared

  cell0.style.width = (cellWidth + 50) + 'px';
  cell0.style.height = cellHeight + 'px';
  cell0.style.backgroundColor = "white";
  cell0.style.position = 'absolute';
  cell0.style.left = (x * cellWidth + 6) + 'px';
  cell0.style.top = '-25px';

  if (view_type == "Species") {
    cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option selected value="Species">Species</option><option value="Protein">Gene Model</option></select>';
  } else {
    cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option value="Species">Species</option><option selected value="Protein">Gene Model</option></select>';
  }

  cell0.style.paddingLeft = "6px";
  zoomedHeatmapEl.appendChild(cell0);

  // -----------------------------
  // Right-hand labels (genomes)
  // -----------------------------
  let save_top = 0; // will become the TOP for the sequence row

  for (let i = 1; i < GN_size; i++) {
    const cell = document.createElement('div');
    x = 50;
    yPosition = i;

    cell.style.width = (cellWidth + 15) + 'px';
    cell.style.height = cellHeight + 'px';
    cell.style.backgroundColor = "white";
    cell.style.position = 'absolute';
    cell.style.left = (x * cellWidth + 6) + 'px';
    cell.style.top  = (cellHeight * (yPosition - 1)) + 'px';

    // temp fix for multiple Mo17 names
    let GN_name_text = String(GN_array[i]);

    if (GM_array[i].startsWith("Zm00014ba")) {
      GN_name_text = "Mo17_v2";
    } else if (GM_array[i].startsWith("Zm00014a")) {
      GN_name_text = "Mo17_v1";
    } else {
      GN_name_text = String(GN_array[i]);
    }

    if (view_type == "Species") {
      cell.innerHTML =
        "<span class='species_option'>" + GN_name_text + "</span>" +
        "<span class='protein_option' style='display: none;'>" + String(GM_array[i]) + "</span>";
    } else {
      // NOTE: your original had a quote bug in class/style; fix it:
      cell.innerHTML =
        "<span class='species_option' style='display: none;'>" + GN_name_text + "</span>" +
        "<span class='protein_option' style='display: inline;'>" + String(GM_array[i]) + "</span>";
    }

    cell.style.paddingLeft = "6px";
    cell.style.color = colorGenome(String(GN_array[i]));
    zoomedHeatmapEl.appendChild(cell);

    // After the last label row, the sequence row should start one row below it
    save_top = (cellHeight * (yPosition)) + 'px';
  }

  // -----------------------------
  // Tooltip setup (unchanged)
  // -----------------------------
  let debounceTimeout2;

  const tooltip_h_pan = document.createElement("div");
  tooltip_h_pan.classList.add("tooltip");
  document.body.appendChild(tooltip_h_pan);

  function showTooltipPan(e) {
    const cell = e.target.closest('.heatmap-cell');
    if (!cell) return;

    const cellData = cell.dataset;
    tooltip_h_pan.style.visibility = "visible";
    tooltip_h_pan.innerHTML = `
      B73 Position: ${cellData.X2}<br>
      Target Position: ${cellData.X3}<br>
      Genome: ${GN_array[cellData.Y]}<br>
      G.M.: ${GM_array[cellData.Y]}<br>
      Substitution: ${cellData.WT} to ${cellData.Sub}<br>
      Score: ${cellData.Score}
    `;
    tooltip_h_pan.style.left = (e.pageX + 10) + "px";
    tooltip_h_pan.style.top  = (e.pageY + 10) + "px";
  }

  zoomedHeatmapEl.addEventListener("mouseover", (e) => {
    if (debounceTimeout2) clearTimeout(debounceTimeout2);
    debounceTimeout2 = setTimeout(() => showTooltipPan(e), 500);
  });

  zoomedHeatmapEl.addEventListener("mouseout", (e) => {
    if (debounceTimeout2) clearTimeout(debounceTimeout2);
    const cell = e.target.closest('.heatmap-cell');
    if (!cell) return;
    tooltip_h_pan.style.visibility = "hidden";
  });

  // -----------------------------
  // ---- Sequence row aligned to HEATMAP-X (no gap columns) ----
  // Build WT lookup from the heatmap data in the current window
  // -----------------------------
  const WT_by_heatX = {};
  data.forEach(d => {
    const hx = +d.X;
    if (hx < start || hx > end) return;

    const wt = String(d.WT ?? "");
    if (!wt || wt === '-') return;

    // keep first WT seen for this column
    if (WT_by_heatX[hx] === undefined) WT_by_heatX[hx] = wt;
  });

  // background strip across the 50-position window
  const seqBg = document.createElement('div');
  seqBg.style.position = 'absolute';
  seqBg.style.left = '0px';
  seqBg.style.top = save_top;
  seqBg.style.width = (50 * cellWidth) + 'px';
  seqBg.style.height = cellHeight + 'px';
  seqBg.style.backgroundColor = 'white';
  zoomedHeatmapEl.appendChild(seqBg);

  // one AA per column, centered
  for (let idx = 0; idx < 50; idx++) {
    const heatX = start + idx;            // HEATMAP coordinate
    const aa = WT_by_heatX[heatX] || "";  // blank if missing

    const letter = document.createElement('div');
    letter.className = 'label_bottom'; // or make a new CSS class
    letter.style.position = 'absolute';
    letter.style.left = (idx * cellWidth) + 'px';
    letter.style.top = save_top;
    letter.style.width = cellWidth + 'px';
    letter.style.height = cellHeight + 'px';
    letter.style.backgroundColor = 'white';
    letter.style.color = 'black';
    letter.style.textAlign = 'center';
    letter.style.lineHeight = cellHeight + 'px';
    letter.style.fontFamily = 'monospace';
    letter.style.fontSize = '14px';
    letter.innerText = aa;

    zoomedHeatmapEl.appendChild(letter);
  }

  // blank cell on the right in the label column
  const seqRight = document.createElement('div');
  seqRight.style.position = 'absolute';
  seqRight.style.left = (50 * cellWidth + 6) + 'px';
  seqRight.style.top = save_top;
  seqRight.style.width = (cellWidth + 15) + 'px';
  seqRight.style.height = cellHeight + 'px';
  seqRight.style.backgroundColor = 'white';
  zoomedHeatmapEl.appendChild(seqRight);
}



    function updateHeatmapZoomPan_OLD(data, start) {
        const end = start + 49;
        const zoomedHeatmapEl = document.getElementById('zoomed-heatmap-pan');
        const cellWidth = window_length / 50;
        const cellHeight = 20;
        const fragment = document.createDocumentFragment();

        zoomedHeatmapEl.innerHTML = ''; // Clear previous zoomed heatmap content

        data.forEach(cellData => {
            const x = +cellData["X"];

            if (x >= start && x <= end) {
                const cell = document.createElement('div');
                cell.className = 'heatmap-cell';
                const xPosition = +cellData["X"];  // + ensures the string is converted to a number
                const yPosition = +cellData["Y"];
                //const score = +cellData["Score"];
                const wild_type = +cellData["WT"];
                const mutation = +cellData["Sub"];
                cell.style.width = (cellWidth + 0.5) + 'px';
                cell.style.height = cellHeight + 'px';
                //cell.style.backgroundColor = colorScalePan(score);
                const wt = String(cellData["WT"] ?? "").trim();
                const sub = String(cellData["Sub"] ?? "").trim();
                const scoreRaw = cellData["Score"];
                const score = (scoreRaw === "" || scoreRaw == null) ? NaN : +scoreRaw;

                const isGap = (wt.startsWith("-") || sub.startsWith("-") || scoreRaw === "" || scoreRaw == null || Number.isNaN(score));

                cell.style.backgroundColor = isGap ? "#e6e6e6" : colorScalePan(score);

                cell.style.position = 'absolute';
                cell.style.left = ((x - start) * cellWidth) + 'px';  // -1 because X is 1-indexed
                cell.style.top = cellHeight * (yPosition - 1) + 'px';  // -1 because Y is 1-indexed
                cell.dataset.X = cellData["X"];
                cell.dataset.X2 = cellData["X2"];
                cell.dataset.Y = cellData["Y"];
                cell.dataset.X3 = cellData["X3"];
                cell.dataset.WT = cellData["WT"];
                cell.dataset.Sub = cellData["Sub"];
                cell.dataset.Score = cellData["Score"];
                fragment.appendChild(cell);
                cell.innerText = cellData["Sub"];
                cell.style.paddingLeft = "6px";
            }
        });

        zoomedHeatmapEl.appendChild(fragment);

        //Drop down menu to toggle between between Species and Protein IDs
        const cell0 = document.createElement('div');
        x = 50
        yPosition = i
        cell0.style.width = (cellWidth + 50) + 'px';
        cell0.style.height = cellHeight + 'px';
        cell0.style.backgroundColor = "white";
        cell0.style.position = 'absolute';
        cell0.style.left = (x * cellWidth + 6) + 'px';
        cell0.style.top =  '-25px';

        if (view_type == "Species") {
            cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option selected value="Species">Species</option><option value="Protein">Gene Model</option></select>';
        } else {
            cell0.innerHTML = '<select id="viewSelector" onchange="switchView()"><option value="Species">Species</option><option selected value="Protein">Gene Model</option></select>';
        }

        cell0.style.paddingLeft = "6px";
        zoomedHeatmapEl.appendChild(cell0);

        let save_top = 0;
        for (let i = 1; i < GN_size; i++) {
            const cell = document.createElement('div');
            x = 50
            yPosition = i
            cell.style.width = (cellWidth + 15) + 'px';
            cell.style.height = cellHeight + 'px';
            cell.style.backgroundColor = "white";
            cell.style.position = 'absolute';
            cell.style.left = (x * cellWidth + 6) + 'px';  // -1 because X is 1-indexed
            cell.style.top = cellHeight * (yPosition - 1) + 'px';  // -1 because Y is 1-indexed

            //temp fix for multiple Mo17 names
            let GN_name_text = String(GN_array[i]);

            if (GM_array[i].startsWith("Zm00014ba"))
            {
                GN_name_text = String("Mo17_v2");
            } else if (GM_array[i].startsWith("Zm00014a"))
            {
                GN_name_text = String("Mo17_v1");
            } else {
                GN_name_text = String(GN_array[i]);
            }

            if (view_type == "Species") {
                cell.innerHTML = "<span class='species_option'>" + GN_name_text + "</span><span class='protein_option' style='display: none;'>" + String(GM_array[i]) + "</span>";
            } else {
                cell.innerHTML = "<span class='species_option' style='display: none;'>" + GN_name_text + "</span><span class='protein_option style='display: inline;'>" + String(GM_array[i]) + "</span>";
            }

            cell.style.paddingLeft = "6px";
            cell.style.color = colorGenome(String(GN_array[i])); // Code that sets the text color based on the Heterotic Group
            zoomedHeatmapEl.appendChild(cell);
            impute = cellHeight * (yPosition ) + 'px';
        }

        //Create the color bar legend
        const colors = [
            '#00429d', '#3860aa', '#587fb3', '#78a0b7', '#9ac0b3', '#c1e19e',
            '#ffff00', '#ffd337', '#fea447', '#f1784d', '#db4c4d', '#bd2147', '#93003a'
        ];

        let debounceTimeout2;

        const tooltip_h_pan = document.createElement("div");
            tooltip_h_pan.classList.add("tooltip");
            document.body.appendChild(tooltip_h_pan);

        function showTooltipPan(e) {
            const cell = e.target.closest('.heatmap-cell');
            if (!cell) return;

            const cellData = cell.dataset;
            tooltip_h_pan.style.visibility = "visible";
            tooltip_h_pan.innerHTML = `
                B73 Position: ${cellData.X2}<br>
                Target Position: ${cellData.X3}<br>
                Genome: ${GN_array[cellData.Y]}<br>
                G.M.: ${GM_array[cellData.Y]}<br>
                Substitution: ${cellData.WT} to ${cellData.Sub}<br>
                Score: ${cellData.Score}
            `;
            tooltip_h_pan.style.left = (e.pageX + 10) + "px";
            tooltip_h_pan.style.top = (e.pageY + 10) + "px";
        }

        zoomedHeatmapEl.addEventListener("mouseover", (e) => {
            // Clear the previous timeout if it exists
            if (debounceTimeout2) {
                clearTimeout(debounceTimeout2);
            }

            // Set a new timeout
            debounceTimeout2 = setTimeout(() => {
                showTooltipPan(e);
            }, 500); // 300ms delay

        });

        zoomedHeatmapEl.addEventListener("mouseout", (e) => {
            // Clear the timeout on mouseout to ensure no tooltips are shown after moving out
            if (debounceTimeout2) {
                clearTimeout(debounceTimeout2);
            }

            const cell = e.target.closest('.heatmap-cell');
            if (!cell) return;

            tooltip_h_pan.style.visibility = "hidden";
        });

        // ---- Sequence row (actual protein sequence) ----
        // save_top is currently the first free row below the last genome label row

        const firstCell = zoomedHeatmapEl.querySelector('.heatmap-cell');
        const originLeft = firstCell ? firstCell.offsetLeft : 0;

        // Use the same visual width as the cells (you use cellWidth + 0.5 for cells)
        const cellW = cellWidth + 0.5;

        // Background strip
        const seqBg = document.createElement('div');
        seqBg.style.position = 'absolute';
        seqBg.style.left = originLeft + 'px';
        seqBg.style.top = save_top;
        seqBg.style.width = (50 * cellW) + 'px';
        seqBg.style.height = cellHeight + 'px';
        seqBg.style.backgroundColor = 'white';
        seqBg.style.zIndex = '5';
        zoomedHeatmapEl.appendChild(seqBg);

        // Letters
        for (let idx = 0; idx < 50; idx++) {
          const alnX = start + idx;
          const aa = (X_to_WT && X_to_WT[String(alnX)]) ? X_to_WT[String(alnX)] : "-";

          const letter = document.createElement('div');
          letter.className = 'label_bottom';
          letter.style.position = 'absolute';      // force absolute, even if CSS says otherwise
          letter.style.left = (originLeft + idx * cellW) + 'px';
          letter.style.top = save_top;
          letter.style.width = cellW + 'px';
          letter.style.height = cellHeight + 'px';
          letter.style.backgroundColor = 'white';
          letter.style.zIndex = '6';
          letter.style.textAlign = 'center';
          letter.style.lineHeight = cellHeight + 'px';
          letter.style.fontFamily = 'monospace';
          letter.style.color = 'black';
          letter.style.paddingLeft = '0px';
          letter.innerText = aa;

          zoomedHeatmapEl.appendChild(letter);
        }


        // Optional: blank cell on the right in the label column so the row looks “complete”
        const seqRight = document.createElement('div');
        seqRight.style.position = 'absolute';
        seqRight.style.left = (50 * cellWidth + 6) + 'px';
        seqRight.style.top = save_top;
        seqRight.style.width = (cellWidth + 15) + 'px';
        seqRight.style.height = cellHeight + 'px';
        seqRight.style.backgroundColor = 'white';
        zoomedHeatmapEl.appendChild(seqRight);


    }

    // Function to create the WT residues on the bottom of the heatmap - this function is no longer needed
    function createZoomedWTLinePan_old(container,data, start) {
        container.innerHTML = "";
        const end = start + 49;

        i = 0
        data.forEach(cellData => {
            const x = +cellData["X"];
            const yPosition = +cellData["Y"]
            const wild_type = cellData["WT"];

        if (x >= start && x <= end && yPosition == "2") {
            const label = document.createElement('div');
            label.className = 'label_bottom';
            label.style.left = scaleZoom(i) + 'px';
            label.innerText = wild_type;
            //label.innerText = x;
            container.appendChild(label);
            i = i + 1;
        }
    });
}

// Function to create the WT (query/B73) residues on the bottom of the zoomed heatmap (alignment-space, gapped)
function createZoomedWTLinePan(container, data, start) {
  container.innerHTML = "";

  const cellWidth = window_length / 50;
  const cellHeight = 20;          // keep consistent with zoom heatmap
  const end = start + 49;

  // One residue per alignment column in the current window
  for (let idx = 0; idx < 50; idx++) {
    const alnX = start + idx;                         // alignment column (1-based)
    const aa = (X_to_WT && X_to_WT[String(alnX)]) ? X_to_WT[String(alnX)] : "-";

    const label = document.createElement("div");
    label.className = "label_bottom";
    label.style.position = "absolute";
    label.style.left = (idx * cellWidth) + "px";
    label.style.top = "0px";
    label.style.width = cellWidth + "px";
    label.style.height = cellHeight + "px";
    label.style.textAlign = "center";
    label.style.lineHeight = cellHeight + "px";
    label.style.fontFamily = "monospace";

    // show '-' for insertion columns relative to query
    label.innerText = aa;

    container.appendChild(label);
  }
}
