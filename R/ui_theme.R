
## shiny application theme
app_theme <- \(){bslib::bs_theme(
  # brand
  brand = FALSE,
  #brand = "_brand.yaml",
  # Controls the default grayscale palette
  bg = "#F8F9FA", fg = "#5200ff",
  # Controls the accent (e.g., hyperlink, button, etc) colors
  primary = "#5200ff", 
  #secondary = "#404040",
  secondary = "#5200ff",
  base_font = c("Grandstander", "sans-serif"),
  code_font = c("Courier", "monospace"),
  heading_font = "'Helvetica Neue', Helvetica, sans-serif",
  # Can also add lower-level customization
  "input-border-color" = "#5200ff"#,
  #"bs-navbar-active-color" = "#c2fe0b",
  #"bs-nav-underline-link-active-color" = "#c2fe0b"
)}



