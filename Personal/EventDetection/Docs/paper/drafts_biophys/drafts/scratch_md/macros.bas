Attribute VB_Name = "NewMacros"
Sub wrap_table()
Attribute wrap_table.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.Macro3"
'
' wrap_table Macro
'
'
    ' Make an index array; only want to align the first N tables
    '
    Selection.WholeStory
    Dim i As Integer
    i = 0
    ' Loop through each table
    '
    Dim t As Table
    For Each t In ActiveDocument.Tables
        i = i + 1
        If (i <= 3) Then
            t.Rows.WrapAroundText = True
            t.Rows.HorizontalPosition = wdTableLeft
            t.Rows.VerticalPosition = InchesToPoints(1)
            t.Rows.RelativeVerticalPosition = wdRelativeVerticalPositionPage
            t.Rows.AllowOverlap = False
            t.Rows.AllowBreakAcrossPages = False
            ' Make first two tables one-column by default
            If (i <= 2) Then
                With t
                    .TopPadding = InchesToPoints(0)
                    .BottomPadding = InchesToPoints(0)
                    .LeftPadding = InchesToPoints(0.02)
                    .RightPadding = InchesToPoints(0.02)
                    .Spacing = 0
                    .AllowPageBreaks = False
                    .AllowAutoFit = False
                End With
                t.PreferredWidthType = wdPreferredWidthPoints
                t.PreferredWidth = InchesToPoints(3.5)
                t.Columns.PreferredWidthType = wdPreferredWidthPoints
                t.Columns.PreferredWidth = InchesToPoints(3.25)
            End If
        End If
    Next
    'Biophysical journal wants...
    ' (1) justified paragraphs
    Selection.ParagraphFormat.Alignment = wdAlignParagraphJustify
    ' (2) Time New Roman Font
    Selection.Font.Name = "Times New Roman"
End Sub
