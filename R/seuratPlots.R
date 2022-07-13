# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


seurat_barplot <- function(object,ident.use,features_list,
                           ident.colors = '',
                           show_percentage_legend=T,
                           percentage_legend_size=2,path_to_save=getwd(),width=10, height=14, ncol=1,plot_name='Bar_plot_percentage_features'){

  require('tidyverse')
  require('dplyr')
  require('ggplot2')
  require('ggpubr')
  require('cowplot')
  require('gridExtra')
  require('RColorBrewer')

  if(unique(ident.colors=='')){
    ident.colors<- c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Pastel2'))
  }

  Idents(object) <- object[[ident.use]]
  idents <- as.data.frame(Idents(object))
  colnames(idents) <- c(ident.use)
  idents$barcode <- row.names(idents)
  count_matrix <- as.data.frame(GetAssayData(object = object, slot = "counts")[paste(features_list,sep=''),])

  data <- as.data.frame(t(count_matrix))
  data$barcode <- rownames(data)
  # print(names(data))
  ident.cluster <-as.data.frame(Idents(object))
  colnames(ident.cluster) <- 'Cell.Type'
  ident.cluster$barcode <- rownames(ident.cluster)
  data.to.plot <- merge(data, ident.cluster,by='barcode')
  data.to.plot$Cell.Type <- as.factor(data.to.plot$Cell.Type)
  #print(table(data.to.plot))
  plot_list <- list()
  for(gene in features_list){
    data_new <- data.to.plot[,c(gene,'Cell.Type')]
    data_new <-  data_new %>%
      filter(get(gene)>0) %>%
      group_by(get('Cell.Type')) %>%
      summarise(cnt=n()) %>%
      mutate(Percentage=cnt/sum(cnt)*100, geneName= gene)

    colnames(data_new) <- c('Cell.type','count','perc','geneName')
    # print(data_new)
    dd <- ggplot(data_new, aes(y=perc, x= Cell.type ,fill =Cell.type)) +
      geom_bar(stat='identity',position = "stack")+
      coord_cartesian(ylim = c(0,105))+
      theme_cowplot()+
      theme(text = element_text(face='bold'),
            axis.title.y =  element_text(size =12,face = 'bold'),
            axis.text.x= element_text(face="bold",
                                      vjust = 1,
                                      hjust = 1,
                                      size = 12,
                                      angle = 45),
            axis.title.x=  element_blank(),
            legend.text = element_text(size = 12,face = 'plain'),
            legend.title = element_text(size =15,face ='bold'),
            # title = element_text(size = 25,face = 'bold'),
            legend.position = "right",
            plot.title = element_text(size =20,face = 'bold') )+
      scale_fill_manual(values = ident.colors)+
      labs(title = gene ,y='Percentage (%)',color="grey60")
    # dd <-  dd +coord_flip()
    if(show_percentage_legend == T){
      dd<- dd+ geom_text(aes(label=round(perc, digits = 2)), vjust=-0.3, size=percentage_legend_size)
    }
    plot_list[[gene]] = dd
  }
  pp <- print(wrap_plots(plot_list)+plot_layout(heights=1, ncol=ncol, guides = 'collect', widths = 2))
  pp <- pp + plot_layout(guides = "collect")
  ggsave(paste(path_to_save,'/',plot_name,'.png',sep=''),plot =pp,width = width, height = height)
  ggsave(paste(path_to_save,'/',plot_name,'.pdf',sep=''),plot =pp,width = width, height = height)
  ggsave(paste(path_to_save,'/',plot_name,'.tiff',sep=''),plot =pp,width = width, height = height)
  print(pp)
  print(paste('save file in : ',path_to_save,'/',plot_name,' (.pdf,.png,.tiff)',sep=''))

}

