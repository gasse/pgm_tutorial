exact.cpt = function(fitted, event, evidence, debug = FALSE){
  
  extract.names = function(call) {
    names = character(0)
    for (x in as.list(call)[-1]) {
      if (is.call(x)) {
        names = c(names, extract.names(x))
      }#THEN
      else if (is.name(x)) {
        names = c(names, as.character(x))
      }#ELSE
    }#FOR
    return(names)
  }#EXTRACT.NAMES
  
  if(is.call(event) | is.logical(event)) {
    nodes = extract.names(event)
  }#THEN
  else {
    nodes = event
  }#ELSE
  
  if(is.call(evidence) | is.logical(evidence)) {
    nodes = union(nodes, extract.names(evidence))
  }#THEN
  else {
    nodes = union(nodes, evidence)
  }#ELSE
  
  nodes = intersect(nodes, names(fitted))
  
  # Joint distribution of the target and conditional nodes knowing their
  # parents, their parent's parents, etc.
  to.check = nodes
  while(length(to.check) > 0) {
    for (node in to.check) {
      parents = fitted[[node]]$parents
      to.check = setdiff(to.check, node)
      to.check = union(to.check, setdiff(parents, nodes))
      nodes = union(nodes, parents)
    }#FOR
  }#WHILE
  
  nbnodes = length(nodes)
  
  cpt.table = expand.grid(lapply(nodes, function(x) {
    dimnames(fitted[[x]]$prob)[[1]]
  }))
  names(cpt.table) = nodes
  
  m = nrow(cpt.table)
  
  # evaluate the expression defining the evidence.
  if (identical(evidence, TRUE))
    r = rep(TRUE, m)
  else
    r = eval(evidence, cpt.table, parent.frame())
  
  # double check that it is a logical vector.
  if (!is.logical(r))
    stop("evidence must evaluate to a logical vector.")
  # double check that it is of the right length.
  if (length(r) != m)
    stop("logical vector for evidence is of length ", length(r), " instead of ", m, ".")
  
  cpt.table = cpt.table[r, , drop=FALSE]
  
  for (node in nodes) {
    cpt = as.data.frame(fitted[[node]]$prob)
    names(cpt)[1] = node # FIX: with 0 parents variable name is missing
    names(cpt)[ncol(cpt)] = paste("p", node, sep=".")
    cpt.table = merge(cpt.table, cpt, by=names(cpt)[-length(cpt)])
  }#FOR
  
  joint.table = cbind(cpt.table[, nodes, drop=FALSE], "p"=1)
  for(i in 1:nbnodes) {
    joint.table[, nbnodes + 1] = joint.table[, nbnodes + 1] * cpt.table[, nbnodes + i]
  }#FOR
  
  return(joint.table)
  
}#EXACT.CPT

exact.dist = function(fitted, event, evidence, debug = FALSE) {
  
  probs = exact.cpt(fitted, event, evidence, debug)
  
  for(col in 1:(ncol(probs)-1)) {
    probs[, col] = factor(probs[, col], exclude=NULL)
    probs = probs[order(probs[, col]), ]
  }#FOR
  
  p.col = ncol(probs)
  cpt = table(probs[, -p.col], dnn=names(probs)[-p.col])
  cpt[1:length(cpt)] = probs[, p.col]
  cpt = prop.table(margin.table(cpt, 1:length(event)))
  
  return(cpt)
  
}#EXACT.DIST
