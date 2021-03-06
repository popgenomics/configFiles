highlight Normal guibg=black guifg=white
set nocompatible
set backspace=indent,eol,start
set ruler
set number
set showcmd
set incsearch
set hlsearch
set hidden 

filetype off

set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

Plugin 'gmarik/Vundle.vim'

" vim en beau
Plugin 'altercation/vim-colors-solarized'
Plugin 'tomasr/molokai'
Plugin 'bling/vim-airline'

" vim en editeur pour programmeur
" Plugin 'scrooloose/nerdtree'
" Plugin 'jistr/vim-nerdtree-tabs'
" Plugin 'klen/python-mode'
 let g:pymode_options_colorcolumn = 0
Plugin 'Valloric/YouCompleteMe'
"let g:ycm_global_ycm_extra_conf = '/home/croux/.vim/bundle/YouCompleteMe/third_party/ycmd/cpp/ycm/.ycm_extra_conf.py'
"let g:ycm_global_ycm_extra_conf = '~/.vim/plugged/youcompleteme/third_party/ycmd/cpp/ycm/.ycm_extra_conf.py'
"let g:ycm_global_ycm_extra_conf = '/home/croux/.vim/bundle/YouCompleteMe/third_party/ycmd/ycmd/tests/clang/testdata/noflags/.ycm_extra_conf.py'
let g:ycm_global_ycm_extra_conf = '/home/croux/.vim/bundle/YouCompleteMe/third_party/ycmd/.ycm_extra_conf.py'

call vundle#end()

" execute pathogen#infect()
syntax on
" filetype plugin indent on

let g:solarized_termcolors=256
set term=xterm-256color
syntax enable
set background=dark
colorscheme default
color desert 

set laststatus=2

let g:airline_detect_paste=1
let g:airline#extensions#tabline#enabled = 1
let g:solarized_termcolors=256
let g:molokai_original = 1

highlight Visual cterm=NONE ctermbg=0 ctermfg=NONE guibg=Grey40

map £ 0



"Copy Paste in the xclipboard
if has("unnamedplus") && has("xterm_clipboard") && hostname() == 'nappa'
set clipboard=unnamedplus
endif


