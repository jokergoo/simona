

#' A shiny app for the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @export
#' @import shiny
#' @importFrom graphics barplot par
#' @examples
#' if(FALSE) {
#'     dag = create_ontology_DAG_from_GO_db()
#'     dag_shiny(dag)
#' }
dag_shiny = function(dag) {

check_pkg("DiagrammeR", bioc = FALSE)
check_pkg("InteractiveComplexHeatmap", bioc = TRUE)

HAS_MCOL = !is.null(mcols(dag))
HAS_NAME_COLUMN = if(HAS_MCOL) "name" %in% colnames(mcols(dag)) else FALSE


ui = fluidPage(
	title = "ontology_DAG browser",
	tags$style(HTML("
body {
	width: 1000px;
}
ul.term_list {
	list-style: circle;
	padding-left:15px;
}

ul.term_list li:hover {
	background-color: #EEEEEE;
}
ul.term_list li {
	background-color: white;
}

#traverse h4 {
	text-align:center;
}
.box {
	float:left;
	border:1px solid #CCCCCC;
	border-radius:4px;
	padding:4px 8px;
	margin-right:20px;
}
.tab-content {
	margin-top:20px;
}
		")),
	h3("ontology_DAG browser"),
	tabsetPanel(
		tabPanel("Summary", 
			div(
				verbatimTextOutput("print_object"),
				plotOutput("summary", width = "1000px", height = "333px")
			)
		),
		tabPanel("Circular Vis",
			plotOutput("circular_plot", width = "900px", height = "700px")
		),
		tabPanel("Traverse", div(
				id = "traverse",
				textInput("term_search", "Term", value = random_terms(dag, 1)),
				actionButton("term_search_submit", "Submit", onclick="$('#term_search_div').show();"),
				div(id = "term_search_div",
					div(
						div(h4("Parent terms"),
							htmlOutput("traverse_parents"),
							class = "box",
							style="width:300px;height:400px;overflow:scroll;"),
						div(h4("Sibling terms"),
							htmlOutput("traverse_siblings"),
							class = "box",
							style="width:300px;height:400px;overflow:scroll;"),
						div(h4("Child terms"),
							htmlOutput("traverse_children"),
							class = "box",
							style="width:300px;height:400px;overflow:scroll;"),
						div(style="clear:both"),
						style = "margin-top:20px;"
					),
					tags$hr(),
					div(
						div(h4("Information"),
							tableOutput("traverse_term_table"), 
							class = "box",
							style="width:300px;"),
						div(h4("Upstream DAG"),
							DiagrammeR::grVizOutput("traverse_term_diagram", width = "100%", height = "550px"), 
							class = "box",
							style="width:600px;height:600px;"),
						div(style="clear:both")
					),
					style = "display:none;"
				)
			)
		),
		tabPanel("Similarity", 
			div(
				div(
					textAreaInput("term_list", "Terms", value = paste(random_terms(dag, min(400, dag_n_terms(dag))), collapse = "\n"), height = 250),
					{	
						if(has_annotation(dag)) {
							all_sim_methods = all_term_sim_methods();
							default_sim_methods = "Sim_Lin_1998"
						} else {
							all_sim_methods = all_term_sim_methods(require_anno = FALSE);
							default_sim_methods = "Sim_WP_1994"
						}
						all_sim_methods = structure(all_sim_methods, names = all_sim_methods)
						radioButtons("term_sim_method", "Select a term similarity method", choices = all_sim_methods,
							selected = default_sim_methods, inline = TRUE)
					},
					tags$script(HTML('
$(".radio-inline").width(150);
$(".radio-inline").first().css("margin-left", 10);
$("#term_sim_method").css("float", "left").css("width", 600).css("margin-left", 20);
$("#term_list").parents("div").css("float", "left");
$(".clear").css("clear", "both");
					'))
				),
				div(class = "clear", style = "clear:both"),
				actionButton("term_list_submit", "Submit", onclick="$('#similarity_div').show();"),
				div(id="similarity_div",
					tabsetPanel(
						tabPanel("Heatmap", 
							InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("similarity_heatmap")
						),
						tabPanel("Circular Vis", 
							plotOutput("circular_plot_highlight", width = "900px", height = "700px")
						)
					),
					style = "margin-top:20px;display:none;"
				)
			),
			tags$script(HTML("
var words = $('#term_list').val().split(/\\s+/);
$('#term_list_submit').text('Submit ' + words.length + ' terms');
document.getElementById('term_list').addEventListener('input', function () {
	var words = $('#term_list').val().split(/\\s+/);
	$('#term_list_submit').text('Submit ' + words.length + ' terms');
})
			"))
		),
		tabPanel("Common ancestors",
			div(
				textInput("term1_ca", "Term 1", value = random_terms(dag, 1)),
				textInput("term2_ca", "Term 2", value = random_terms(dag, 1)),
				actionButton("term_ca_submit", "Submit"),
				DiagrammeR::grVizOutput("term_ca_diagram", width = "600px", height = "600px")
			)
		)
	),
	tags$script(HTML('
var target = document.getElementById("traverse_term_diagram");

var observer = new MutationObserver(function (mutations) {
    mutations.forEach(function (mutation) {
        $("#traverse_term_diagram svg g.node").hover(function() {
			$( this ).css("cursor", "pointer").click(function() {
				var term_id = $( this ).find("title").text();
				Shiny.onInputChange("traverse_term", term_id);
			})
		})
    });
});

var config = {
    attributes: true,
    childList: true,
    characterData: true
};

observer.observe(target, config);
	')),

	hr(),
	HTML(qq('<p style="clear:both">Generated by <a href="https://bioconductor.org/packages/simona/" target="_blank">simona @{installed.packages()["simona", "Version"]}</a>.</p>')),
)

server = function(input, output, session) {
	output$print_object = renderPrint({
		print(dag)
	})

	output$circular_plot = renderPlot({
		showNotification("Generating circular vis...", duration = 5)
		dag_circular_viz(dag, partition_by_size = round(dag_n_terms(dag)/5))
	})

	output$summary = renderPlot({
		np = n_parents(dag)
		nc = n_children(dag)
		dp = dag_depth(dag)

		par(mfrow = c(1, 3))
		barplot(table(np), xlab = "Numbers of parents", ylab = "Counts")
		barplot(table(nc), xlab = "Numbers of children", ylab = "Counts")
		barplot(table(dp), xlab = "Depths", ylab = "Counts")
	})

	as_li = function(term_id, x = term_id, text = "") {
		if(length(term_id) == 0) {
			html = text
		} else {
			html = qq("<li><a title='@{term_id}' style='color:#337ab7;text-decoration:underline;cursor:pointer;' onclick='Shiny.onInputChange(\"traverse_term\", \"@{term_id}\")'>@{x}</a></li>\n")
			html = paste0("<ul class='term_list'>\n", html, "\n</ul>")
		}
		HTML(html)
	}

	update_travese = function(dag, term) {

		if(!dag_has_terms(dag, term)) {
			throw_error(qq("Term '@{term}' does not exist in the ontology."))
			return(NULL)
		}

		updateTextInput(inputId = "term_search", value = term)

		parents = dag_parents(dag, term)
		parents = setdiff(parents, term)

		siblings = dag_siblings(dag, term)
		siblings = c(term, siblings)

		children = dag_children(dag, term)
		children = setdiff(children, term)

		if(HAS_NAME_COLUMN) {
			mc = mcols(dag)
			parent_names = ifelse(is.na(mc[parents, "name"]), parents, mc[parents, "name"])
			sibling_names = ifelse(is.na(mc[siblings, "name"]), siblings, mc[siblings, "name"])
			child_names = ifelse(is.na(mc[children, "name"]), children, mc[children, "name"])
		} else {
			parent_names = parents
			sibling_names = siblings
			child_names = children
		}

		output$traverse_parents = renderUI({
			as_li(parents, parent_names, "Already root, no parent.")
		})
		output$traverse_siblings = renderUI({
			sibling_names[1] = paste0("<b style='color:red'>", sibling_names[1], "</b>")
			as_li(siblings, sibling_names)
		})
		output$traverse_children = renderUI({
			as_li(children, child_names, "Already leaf, no child.")
		})

		output$traverse_term_table = renderTable({
			tb = mcols(dag)
			if(is.null(tb)) {
				tb = data.frame(id = term)
			} else {
				tb = tb[term, ,drop = FALSE]
			}

			tb = t(tb)
			tb[is.na(tb)] = ""
			tb
		}, rownames = TRUE, colnames = FALSE)

		output$traverse_term_diagram = DiagrammeR::renderGrViz({
			dag_graphviz(dag[, term], node_param = list(fillcolor = structure(names = term, "red"), style = "filled"), width = 600, height = 600)
		})
	}

	observeEvent(input$term_search_submit, {
		term = input$term_search
		showNotification(qq("Search term '@{term}'..."), duration = 2)
		update_travese(dag, term)
	})

	observeEvent(input$traverse_term, {
		term = input$traverse_term
		showNotification(qq("Search term '@{term}'..."), duration = 2)
		update_travese(dag, term)
	})

	observeEvent(input$term_list_submit, {
		showNotification("Generating heatmap...", duration = 5)
		
		terms = strsplit(input$term_list, "\\s+")[[1]]
		terms = terms[!grepl("^\\s*$", terms)]

		terms = terms[dag_has_terms(dag, terms)]

		if(length(terms) < 2) {
			throw_error(qq("Terms do not exist in the ontology."))
		}

		method = input$term_sim_method
		sim = term_sim(dag, terms, method = method)
		ht = Heatmap(sim, name = "Similarity", show_row_names = FALSE, show_column_names = FALSE,
			show_row_dend = FALSE, show_column_dend = FALSE, column_title = qq("@{method} on @{nrow(sim)} terms"))
		InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, ht, "similarity_heatmap")

		output$circular_plot_highlight = renderPlot({
			showNotification("Generating circular vis...", duration = 5)
			dag_circular_viz(dag, highlight = terms, partition_by_size = round(dag_n_terms(dag)/5))
		})
		
	})

	observeEvent(input$term_ca_submit, {
		if(input$term1_ca != "" && input$term2_ca != "") {
			term1 = input$term1_ca
			term2 = input$term2_ca

			if(!dag_has_terms(dag, term1)) {
				throw_error(qq("Term '@{term1}' does not exist in the ontology."))
				return(NULL)
			}

			if(!dag_has_terms(dag, term2)) {
				throw_error(qq("Term '@{term2}' does not exist in the ontology."))
				return(NULL)
			}

			showNotification(qq("Search common ancestors of '@{term1}' and '@{term2}'..."), duration = 5)
		

			a1 = dag_ancestors(dag, term1, include_self = TRUE)
			a2 = dag_ancestors(dag, term2, include_self = TRUE)
			aa = union(a1, a2)
			ca = intersect(a1, a2)
				
			fillcolor = structure(names = ca, rep("lightblue", length(ca)))
			fillcolor[term1] = "red"
			fillcolor[term2] = "red"
			node_param = list(fillcolor = fillcolor, style = "filled")
			output$term_ca_diagram = DiagrammeR::renderGrViz({
				dag_graphviz(dag_filter(dag, terms = aa), node_param = node_param, width = 600, height = 600)
			})
		}
	})
}

throw_error = function(message) {
	showModal(modalDialog(
        title = "error",
        p(message, style = "color:red;"),
        easyClose = FALSE,
        size = "m"
    ))
}


print(shinyApp(ui, server))



}


